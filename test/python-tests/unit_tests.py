import numpy as np
import zonoopt as zono
import scipy.sparse as sp
from pathlib import Path

"""zonoLAB used to generate these unit tests"""

# globals: unit test folder
unit_test_folder = Path(__file__).parent

# unit tests
def test_vrep_2_hz():

    # folder where unit test data resides
    test_folder = unit_test_folder / 'vrep_2_hybzono'

    # build hybzono from vrep in zonocpp
    V_polys = []
    V_polys.append(np.array([[5.566, 5.896],
                             [4.044, 5.498],
                             [5.32, 3.909],
                             [5.599, 4.082]]))
    V_polys.append(np.array([[0.049, 6.05],
                             [-0.248, 3.881],
                             [0.617, 3.981]]))
    V_polys.append(np.array([[5.481, 0.911],
                             [4.937, 1.183],
                             [5.199, -1.001]]))
    V_polys.append(np.array([[3.447, 3.207],
                             [2.853, 3.552],
                             [3.341, 1.914],
                             [3.656, 2.397]]))
        
    Z = zono.vrep_2_hybzono(V_polys)

    if Z.is_0_1_form():
        Z.convert_form()

    # expected result
    Gc_expected = np.loadtxt(test_folder / 'Gc.txt', delimiter=' ')
    Gb_expected = np.loadtxt(test_folder / 'Gb.txt', delimiter=' ')
    c_expected = np.loadtxt(test_folder / 'c.txt', delimiter=' ')
    Ac_expected = np.loadtxt(test_folder / 'Ac.txt', delimiter=' ')
    Ab_expected = np.loadtxt(test_folder / 'Ab.txt', delimiter=' ')
    b_expected = np.loadtxt(test_folder / 'b.txt', delimiter=' ')

    # correct equality constraints for normalization
    for i in range(len(b_expected)):
        if np.abs(b_expected[i]-Z.get_b()[i]) > 1e-3:
            Ac_expected[i,:] *= Z.get_b()[i]/b_expected[i]
            Ab_expected[i,:] *= Z.get_b()[i]/b_expected[i]
            b_expected[i] = Z.get_b()[i]

    # compare results
    assert np.allclose(Z.get_Gc().toarray(), Gc_expected)
    assert np.allclose(Z.get_Gb().toarray(), Gb_expected)
    assert np.allclose(Z.get_c(), c_expected)
    assert np.allclose(Z.get_Ac().toarray(), Ac_expected)
    assert np.allclose(Z.get_Ab().toarray(), Ab_expected)
    assert np.allclose(Z.get_b(), b_expected)
    print('Passed: V-rep to Hybzono')

def test_minkowski_sum():

    # folder where unit test data resides
    test_folder = unit_test_folder / 'minkowski_sum'

    # build hybzono from vrep in zonocpp
    V_polys = []
    V_polys.append(np.array([[5.566, 5.896],
                             [4.044, 5.498],
                             [5.32, 3.909],
                             [5.599, 4.082]]))
    V_polys.append(np.array([[0.049, 6.05],
                             [-0.248, 3.881],
                             [0.617, 3.981]]))
    V_polys.append(np.array([[5.481, 0.911],
                             [4.937, 1.183],
                             [5.199, -1.001]]))
    V_polys.append(np.array([[3.447, 3.207],
                             [2.853, 3.552],
                             [3.341, 1.914],
                             [3.656, 2.397]]))
        
    Z1 = zono.vrep_2_hybzono(V_polys)

    # zonotope
    G = 0.5*np.array([[np.sqrt(3), 1, np.sqrt(3)],
                        [0.5, 0, -0.5]])
    c = np.array([-2.0, 1.0])
    Z2 = zono.Zono(G,c)

    # minkowski sum
    Z = zono.minkowski_sum(Z1, Z2)
    if Z.is_0_1_form():
        Z.convert_form()

    # expected result
    Gc_expected = np.loadtxt(test_folder / 'Gc.txt', delimiter=' ')
    Gb_expected = np.loadtxt(test_folder / 'Gb.txt', delimiter=' ')
    c_expected = np.loadtxt(test_folder / 'c.txt', delimiter=' ')
    Ac_expected = np.loadtxt(test_folder / 'Ac.txt', delimiter=' ')
    Ab_expected = np.loadtxt(test_folder / 'Ab.txt', delimiter=' ')
    b_expected = np.loadtxt(test_folder / 'b.txt', delimiter=' ')

    # correct equality constraints for normalization
    for i in range(len(b_expected)):
        if np.abs(b_expected[i]-Z.get_b()[i]) > 1e-3:
            Ac_expected[i,:] *= Z.get_b()[i]/b_expected[i]
            Ab_expected[i,:] *= Z.get_b()[i]/b_expected[i]
            b_expected[i] = Z.get_b()[i]

    # compare results
    assert np.allclose(Z.get_Gc().toarray(), Gc_expected)
    assert np.allclose(Z.get_Gb().toarray(), Gb_expected)
    assert np.allclose(Z.get_c(), c_expected)
    assert np.allclose(Z.get_Ac().toarray(), Ac_expected)
    assert np.allclose(Z.get_Ab().toarray(), Ab_expected)
    assert np.allclose(Z.get_b(), b_expected)
    print('Passed: Minkowski Sum')
    
def test_intersection():

    # folder where unit test data resides
    test_folder = unit_test_folder / 'intersection'

    # build hybzono from vrep in zonocpp
    V_polys = []
    V_polys.append(np.array([[5.566, 5.896],
                             [4.044, 5.498],
                             [5.32, 3.909],
                             [5.599, 4.082]]))
    V_polys.append(np.array([[0.049, 6.05],
                             [-0.248, 3.881],
                             [0.617, 3.981]]))
    V_polys.append(np.array([[5.481, 0.911],
                             [4.937, 1.183],
                             [5.199, -1.001]]))
    V_polys.append(np.array([[3.447, 3.207],
                             [2.853, 3.552],
                             [3.341, 1.914],
                             [3.656, 2.397]]))
        
    Z1 = zono.vrep_2_hybzono(V_polys)

    # zonotope
    G = 0.5*np.array([[np.sqrt(3), 1, np.sqrt(3)],
                        [0.5, 0, -0.5]])
    c = np.array([-2.0, 1.0])
    Z2 = zono.Zono(sp.csc_matrix(G), c)

    # minkowski sum
    Z3 = zono.minkowski_sum(Z1, Z2)

    # conzono
    G = np.array([[3.0, 0.0, 0.0], 
                  [0.0, 3.0, 0.0]])
    c = np.array([-0.5, 4.5])
    A = np.ones((1, 3))
    b = np.array([1.0])
    Z4 = zono.ConZono(sp.csc_matrix(G), c, sp.csc_matrix(A), b)

    # intersection
    Z = zono.intersection(Z3, Z4)
    if Z.is_0_1_form():
        Z.convert_form()

    # expected result
    Gc_expected = np.loadtxt(test_folder / 'Gc.txt', delimiter=' ')
    Gb_expected = np.loadtxt(test_folder / 'Gb.txt', delimiter=' ')
    c_expected = np.loadtxt(test_folder / 'c.txt', delimiter=' ')
    Ac_expected = np.loadtxt(test_folder / 'Ac.txt', delimiter=' ')
    Ab_expected = np.loadtxt(test_folder / 'Ab.txt', delimiter=' ')
    b_expected = np.loadtxt(test_folder / 'b.txt', delimiter=' ')

    # correct equality constraints for normalization
    for i in range(len(b_expected)):
        if np.abs(b_expected[i]-Z.get_b()[i]) > 1e-3:
            Ac_expected[i,:] *= Z.get_b()[i]/b_expected[i]
            Ab_expected[i,:] *= Z.get_b()[i]/b_expected[i]
            b_expected[i] = Z.get_b()[i]

    # compare results
    assert np.allclose(Z.get_Gc().toarray(), Gc_expected)
    assert np.allclose(Z.get_Gb().toarray(), Gb_expected)
    assert np.allclose(Z.get_c(), c_expected)
    assert np.allclose(Z.get_Ac().toarray(), Ac_expected)
    assert np.allclose(Z.get_Ab().toarray(), Ab_expected)
    assert np.allclose(Z.get_b(), b_expected)
    print('Passed: Intersection')

def test_is_empty():

    # folder where unit test data resides
    test_folder = unit_test_folder / 'conzono_feasibility'

    # load in feasible conzono
    G = np.loadtxt(test_folder / 'f_G.txt', delimiter=' ')
    c = np.loadtxt(test_folder / 'f_c.txt', delimiter=' ')
    A = np.loadtxt(test_folder / 'f_A.txt', delimiter=' ')
    b = np.loadtxt(test_folder / 'f_b.txt', delimiter=' ')

    Zf = zono.ConZono(sp.csc_matrix(G), c, sp.csc_matrix(A), b)

    # load in infeasible conzono
    G = np.loadtxt(test_folder / 'i_G.txt', delimiter=' ')
    c = np.loadtxt(test_folder / 'i_c.txt', delimiter=' ')
    A = np.loadtxt(test_folder / 'i_A.txt', delimiter=' ')
    b = np.loadtxt(test_folder / 'i_b.txt', delimiter=' ')

    Zi = zono.ConZono(sp.csc_matrix(G), c, sp.csc_matrix(A), b)
    
    # check if empty
    assert not Zf.is_empty()
    assert Zi.is_empty()
    print('Passed: Is Empty')

def test_support():

    # folder where unit test data resides
    test_folder = unit_test_folder / 'support'

    # load in conzono
    G = np.loadtxt(test_folder / 'G.txt', delimiter=' ')
    c = np.loadtxt(test_folder / 'c.txt', delimiter=' ')
    A = np.loadtxt(test_folder / 'A.txt', delimiter=' ')
    b = np.loadtxt(test_folder / 'b.txt', delimiter=' ')

    Z = zono.ConZono(sp.csc_matrix(G), c, sp.csc_matrix(A), b)

    # load direction and expected support function
    d = np.loadtxt(test_folder / 'd.txt', delimiter=' ')
    s_expected = np.loadtxt(test_folder / 'sup.txt', delimiter=' ')

    # compute support function
    s = Z.support(d)
    
    # compare results
    tol = 5e-2 # tolerance on success
    assert np.abs(s-s_expected)/np.abs(s_expected) < tol
    print('Passed: Support Function')

def test_point_contain():

    # folder where the data resides
    test_folder = unit_test_folder / 'point_containment'

    # load in conzono
    G = np.loadtxt(test_folder / 'G.txt', delimiter=' ')
    c = np.loadtxt(test_folder / 'c.txt', delimiter=' ')
    A = np.loadtxt(test_folder / 'A.txt', delimiter=' ')
    b = np.loadtxt(test_folder / 'b.txt', delimiter=' ')

    Z = zono.ConZono(sp.csc_matrix(G), c, sp.csc_matrix(A), b)

    # load point in set
    x_c = np.loadtxt(test_folder / 'x_c.txt', delimiter=' ')

    # load point not in set
    x_n = np.loadtxt(test_folder / 'x_n.txt', delimiter=' ')

    # check correct classification of containment
    assert Z.contains_point(x_c)
    assert not Z.contains_point(x_n)
    print('Passed: Point Containment')

def test_get_leaves():

    # make random conzonos
    np.random.seed(0)

    n_CZs = 20
    n = 10
    nV = 2*n

    CZs = []
    for i in range(n_CZs):
        V = np.random.random((nV, n))
        CZs.append(zono.vrep_2_conzono(V))

    # take union
    U = zono.union_of_many(CZs)

    # minkowski sum
    Z = zono.minkowski_sum(U, U)

    # get number of leaves
    settings = zono.OptSettings()
    settings.n_threads_bnb = 1
    leaves = Z.get_leaves(settings=settings)

    # check number of leaves is correct
    assert len(leaves) == n_CZs**2
    print('Passed: Get Leaves')

# run the unit tests
test_vrep_2_hz()
test_minkowski_sum()
test_intersection()
test_is_empty()
test_support()
test_point_contain()
test_get_leaves()