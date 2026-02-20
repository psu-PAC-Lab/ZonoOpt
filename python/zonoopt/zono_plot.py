import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
from scipy.spatial import ConvexHull
from scipy.optimize import linprog
from scipy.linalg import null_space
import time
import warnings

from ._core import *


def find_vertex(Z, d):
    """Get vertex of Z nearest to direction d"""
    
    # maximize dot product
    c = -Z.get_G().transpose().dot(d)
    if Z.is_0_1_form():
        bounds = [(0, 1) for i in range(Z.get_nG())]
    else:
        bounds = [(-1, 1) for i in range(Z.get_nG())]

    if Z.is_zono():
        res = linprog(c, bounds=bounds)
    elif Z.is_conzono():
        res = linprog(c, A_eq=Z.get_A(), b_eq=Z.get_b(), bounds=bounds)
    else:
        raise ValueError('find_vertex unsupported data type')

    if res.success:
        return Z.get_G()*res.x + Z.get_c()
    else:
        return None

def get_conzono_vertices(Z, t_max=60.0):
    """Get vertices of Z"""

    # init time
    t0 = time.time()

    # make sure Z is not empty
    if Z.is_empty():
        warnings.warn('Z is empty, returning empty list of vertices.')
        return []

    # search for vertices along perpendicular directions to get initial simplex
    verts = []
    D = [np.array([1 if i==j else 0 for j in range(Z.get_n())]) for i in range(Z.get_n())] # init directions as standard basis
    B = [] # init basis
    for _ in range(Z.get_n()):

        # support in positive direction
        d = D.pop()
        vd_pos = find_vertex(Z, d)    
        
        # support in negative direction
        vd_neg = find_vertex(Z, -d)

        # make sure feasible
        if vd_pos is None or vd_neg is None: # infeasible, not detected during get_leaves
            return []
        
        # check if vertices are new and whether direction is thin
        is_vd_pos_new = not any(np.allclose(vd_pos, v) for v in verts)
        is_vd_neg_new = not any(np.allclose(vd_neg, v) for v in verts)
        is_thin = np.allclose(vd_pos, vd_neg)

        # add new vertices
        if is_vd_pos_new:
            verts.append(vd_pos)
        if is_vd_neg_new and not is_thin:
            verts.append(vd_neg)

        # update direction
        if not is_thin:
            B.append(vd_pos - vd_neg)
            N = null_space(np.array(B)).transpose()
            D = [N[j,:].flatten() for j in range(N.shape[0])]

    # return if set is not full-dimensional
    if len(verts) < Z.get_n()+1:
        return np.array(verts)
    
    # search for additional vertices along the directions of the facet normals
    converged = False
    while not converged and ((time.time()-t0) < t_max):

        # compute convex hull and centroid
        verts_np_arr = np.array(verts)
        hull = ConvexHull(verts_np_arr)
        centroid = np.mean(verts_np_arr, axis=0)

        # get facet normals
        normals = []
        for simplex in hull.simplices:
            
            # get vertices of facet. each row is a vertex
            V = verts_np_arr[simplex]
            
            # get normal
            Vn = V[-1,:] # last element
            A = V[:-1,:] - Vn # subtract last element from each row
            _, _, Vt = np.linalg.svd(A) # singular value decomp to get null space
            n = Vt[-1,:] # last row of Vt is the null space

            # ensure outward normal
            if np.dot(n, Vn - centroid) < 0:
                n = -n

            normals.append(n)

        # search facet normals for additional vertices
        n_new_verts = 0 # init
        for n in normals:

            # get vertex
            vd = find_vertex(Z, n)

            # check if vertex is new
            if not any(np.allclose(vd, v) for v in verts):
                verts.append(vd)
                n_new_verts += 1

        # check for convergence
        if n_new_verts == 0:
            converged = True

    # throw warning if time limit was reached
    if (time.time()-t0) > t_max:
        warnings.warn('get_vertices time limit reached, terminating early.')

    V = np.array(verts)
    hull = ConvexHull(V)
    V = V[hull.vertices,:]

    return V

def get_vertices(Z, t_max=60.0):
    """
    Get vertices of zonotopic set using scipy linprog.
    
    Args:
        Z (HybZono): Zonotopic set.
        t_max (float, optional): Maximum time to spend on finding vertices. Defaults to 60.0 seconds.
    
    Returns:
        numpy.ndarray: Vertices of the zonotopic set. If Z is a point, returns its coordinates.
    """
    
    if Z.is_empty_set():
        return None
    elif Z.is_point():
        return Z.get_c().reshape(1,-1)
    elif Z.is_zono() or Z.is_conzono():
        return get_conzono_vertices(Z, t_max=t_max)
    elif Z.is_hybzono():
        raise ValueError('get_vertices not implemented for HybZono')

def plot(Z, ax=None, settings=OptSettings(), t_max=60.0, **kwargs):
    """
    Plots zonotopic set using matplotlib.

    Args:
        Z (HybZono): zonotopic set to be plotted
        ax (matplotlib.axes.Axes, optional): Axes to plot on. If None, current axes are used.
        settings (OptSettings, optional): Settings for the optimization. Defaults to OptSettings().
        t_max (float, optional): Maximum time to spend on finding vertices. Defaults to 60.0 seconds.
        **kwargs: Additional keyword arguments passed to the plotting function (e.g., color, alpha).

    Returns:
        list: List of matplotlib objects representing the plotted zonotope.
    """

    if Z.get_n() < 2 or Z.get_n() > 3:
        raise ValueError("Plot only implemented in 2D or 3D")
    
    # hybzono -> get leaves
    if Z.is_hybzono():
        leaves = Z.get_leaves(settings=settings)
        if len(leaves) > 0:
            time_per_leaf = t_max / len(leaves)
        else:
            warnings.warn('No leaves found in HybZono, returning empty plot')
            return []
        objs = []
        for leaf in leaves:
            objs.append(plot(leaf, ax=ax, t_max=time_per_leaf, **kwargs)[0])
        return objs

    V = get_vertices(Z, t_max=t_max)

    # 2D
    if Z.get_n() == 2:
        
        # get axes
        if ax is None:
            ax = plt.gca()

        # plot
        if V is None or len(V) == 0:
            warnings.warn("No vertices found, returning empty plot")
            return ax.plot([], [])
        elif Z.is_point() or len(V) < Z.get_n()+1:
            try:
                return ax.plot(V[:,0], V[:,1], **kwargs)
            except Exception as e:
                print(V)
                warnings.warn(f"Error plotting point / line: {e}")
                return ax.plot([], [])
        else:
            return ax.fill(V[:,0], V[:,1], **kwargs)

    else: # 3D

        # get axes
        if ax is None:
            raise ValueError("3D plotting requires an Axes3D object")
        
        # plot
        if V is None or len(V) == 0:
            warnings.warn("No vertices found, returning empty plot")
            obj = ax.scatter([], [], [])
        elif Z.is_point():
            obj = ax.scatter(V[0,0], V[0,1], V[0,2], **kwargs)
        else:
            hull = ConvexHull(V)
            obj = ax.add_collection3d(Poly3DCollection([[V[vertex] for vertex in face] for face in hull.simplices], **kwargs))
            xmin, xmax = ax.get_xlim()
            ymin, ymax = ax.get_ylim()
            zmin, zmax = ax.get_zlim()
            ax.auto_scale_xyz(np.hstack([V[:,0].flatten(), [xmin, xmax]]),
                            np.hstack([V[:,1].flatten(), [ymin, ymax]]),
                            np.hstack([V[:,2].flatten(), [zmin, zmax]]))

        # adjust scaling
        return obj

        
    