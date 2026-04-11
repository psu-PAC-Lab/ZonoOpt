import numpy as np
import time
import warnings

from ._core import *

def get_vertices(Z, t_max=60.0):
    """
    Get vertices of zonotopic set using scipy linprog.
    
    Args:
        Z (HybZono): Zonotopic set.
        t_max (float, optional): Maximum time to spend on finding vertices. Defaults to 60.0 seconds.
    
    Returns:
        numpy.ndarray: Vertices of the zonotopic set. If Z is a point, returns its coordinates.
    """

    # local import to reduce zonoopt import time if plotting is not used
    from scipy.spatial import ConvexHull
    from scipy.optimize import linprog
    from scipy.linalg import null_space

    def _find_vertex(Z, d):
        """Get vertex of Z nearest to direction d"""
        
        # maximize dot product
        c = -Z.get_G().transpose().dot(d)
        if Z.is_0_1_form():
            bounds = [(0, 1) for _ in range(Z.get_nG())]
        else:
            bounds = [(-1, 1) for _ in range(Z.get_nG())]

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

    def _get_conzono_vertices(Z, t_max=60.0):
        """Get vertices of Z"""

        # init time
        t0 = time.time()

        # search for vertices along perpendicular directions to get initial simplex
        verts = np.zeros((0, Z.get_n()))
        D = [np.array([1 if i==j else 0 for j in range(Z.get_n())]) for i in range(Z.get_n())] # init directions as standard basis
        B = [] # init basis
        vc = np.zeros(Z.get_n()) # init vertex candidate as origin
        for _ in range(Z.get_n()):

            # support in positive direction
            d = D.pop()
            vd_pos = _find_vertex(Z, d)    
            
            # support in negative direction
            vd_neg = _find_vertex(Z, -d)

            # make sure feasible
            if vd_pos is None or vd_neg is None: # infeasible, not detected during get_leaves
                return []
            
            # check if vertices are new and whether direction is thin
            is_vd_pos_new = not any(np.allclose(vd_pos, v) for v in verts)
            is_vd_neg_new = not any(np.allclose(vd_neg, v) for v in verts)
            is_thin = np.abs((vd_pos-vc).dot(d)) < 1e-6 and np.abs((vd_neg-vc).dot(-d)) < 1e-6 or np.allclose(vd_pos, vd_neg)
            
            # add new vertices
            if not is_thin:

                # add new vertices
                if is_vd_pos_new:
                    verts = np.vstack([verts, vd_pos])
                if is_vd_neg_new:
                    verts = np.vstack([verts, vd_neg])

                # update direction
                B.append(vd_pos - vd_neg)
                N = null_space(np.array(B)).transpose()
                D = [N[j,:].flatten() for j in range(N.shape[0])]
                vc = (vd_pos + vd_neg) / 2.

        # return if set is not full-dimensional
        if len(verts) < Z.get_n()+1:
            return verts
        
        # search for additional vertices along the directions of the facet normals
        converged = False
        normals = []
        while not converged and ((time.time()-t0) < t_max):

            # compute convex hull and centroid
            verts_np_arr = np.array(verts)
            try:
                hull = ConvexHull(verts_np_arr)
            except:
                warnings.warn('ConvexHull failed, returning current vertices')
                return np.array(verts)
            centroid = np.mean(verts_np_arr, axis=0)

            # get facet normals
            new_normals = []
            for simplex in hull.simplices:
                
                # get vertices of facet. each row is a vertex
                V = verts_np_arr[simplex]
                
                # get normal
                Vn = V[-1,:] # last element
                A = V[:-1,:] - Vn # subtract last element from each row
                N = null_space(A).transpose()
                n = N[0,:]

                # ensure outward normal
                if np.dot(n, Vn - centroid) < 0:
                    n = -n

                if not any(np.allclose(n, existing_n) for existing_n in normals):
                    new_normals.append(n)

            # search facet normals for additional vertices
            n_new_verts = 0 # init
            for n in new_normals:

                # get vertex
                vd = _find_vertex(Z, n)
                if vd is None:
                    warnings.warn('find_vertex failed, skipping direction')
                    continue

                # check if vertex is new
                if not any(np.allclose(vd, v) for v in verts):
                    verts = np.vstack([verts, vd])
                    n_new_verts += 1

            # already-checked normal directions
            normals.extend(new_normals)

            # check for convergence
            converged = n_new_verts == 0


        # throw warning if time limit was reached
        if (time.time()-t0) > t_max:
            warnings.warn('get_vertices time limit reached, terminating early.')

        V = np.array(verts)
        hull = ConvexHull(V)
        V = V[hull.vertices,:]

        return V
    
    # get vertices based on type
    if Z.is_empty_set():
        return np.zeros((0, Z.get_n()))
    elif Z.is_point():
        return Z.get_c().reshape(1,-1)
    elif Z.is_zono() or Z.is_conzono():
        return _get_conzono_vertices(Z, t_max=t_max)
    elif Z.is_hybzono():
        t0 = time.time()
        settings = OptSettings()
        settings.t_max = t_max
        settings.verbose = True
        sol = OptSolution()
        Z_leaves = Z.get_leaves(settings=settings, solution=sol)
        if not sol.converged and not sol.infeasible:
            warnings.warn('get_leaves returned before convergence, get_vertices may be incomplete.')
        dt = time.time() - t0
        V = np.zeros((0, Z.get_n()))
        for leaf in Z_leaves:
            V_leaf = get_vertices(leaf, t_max=t_max-dt)
            if V_leaf is not None:
                V = np.vstack((V, V_leaf))
            dt = time.time() - t0
        return V
    else:
        raise ValueError('get_vertices unsupported data type')

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

    # local import to reduce zonoopt import time if plotting is not used
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from scipy.spatial import ConvexHull
    from tqdm import tqdm

    if Z.get_n() < 2 or Z.get_n() > 3:
        raise ValueError("Plot only implemented in 2D or 3D")
    
    # hybzono -> get leaves
    if Z.is_hybzono():
        t0 = time.time()
        
        sol = OptSolution()
        leaves = Z.get_leaves(settings=settings, solution=sol)

        if not sol.converged and not sol.infeasible:
            warnings.warn('get_leaves returned before convergence, plot may be incomplete.')

        if len(leaves) == 0:
            warnings.warn('No leaves found in HybZono, returning empty plot')
            return []
        
        objs = []
        pbar = tqdm(leaves)
        for leaf in pbar:
            pbar.set_description('Plotting HybZono leaves')
            t = time.time() - t0
            if t > t_max:
                warnings.warn('Plotting time limit reached, terminating early.')
                break
            obj = plot(leaf, ax=ax, t_max=t_max-t, settings=settings, **kwargs)
            if Z.get_n() == 2:
                objs.append(obj[0])
            else:
                objs.append(obj)
        return objs

    V = get_vertices(Z, t_max=t_max)

    # 2D
    if Z.get_n() == 2:
        
        # get axes
        if ax is None:
            ax = plt.gca()

        # plot
        if V is None or len(V) == 0:
            return ax.plot([], [])
        elif len(V) <= 2: # line or point
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
            obj = ax.scatter([], [], [])
        elif len(V) == 1: # point
            obj = ax.scatter(V[:,0], V[:,1], V[:,2], **kwargs)
        elif len(V) == 2: # line
            obj = ax.plot(V[:,0], V[:,1], V[:,2], **kwargs)
        elif len(V) == 3: # plane
            obj = ax.add_collection3d(Poly3DCollection([V], **kwargs))
        else:
            try:
                hull = ConvexHull(V)
                obj = ax.add_collection3d(Poly3DCollection([[V[vertex] for vertex in face] for face in hull.simplices], **kwargs))
            except Exception as e:
                print(V)
                warnings.warn(f"Error plotting 3D zonotope: {e}")
                obj = ax.scatter([], [], [])
        ax.autoscale_view()

        # adjust scaling
        return obj

        
    