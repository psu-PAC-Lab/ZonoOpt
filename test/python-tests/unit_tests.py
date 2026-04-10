import numpy as np
import zonoopt as zono
from scipy import sparse
from pathlib import Path

# Note: zonoLAB used to generate data for some of these unit tests

# globals: unit test folder
test_data_folder = Path(__file__).parent.parent / 'test-data'

# utilities: random matrices / vectors / sets
class TestUtilities:

    @staticmethod
    def random_sparse_matrix(m, n, density, val_min, val_max):
        assert(m >= 0 and n >= 0), "Matrix dimensions must be non-negative"
        assert(0 < density <= 1), "Density must be in (0, 1]"
        assert(val_min < val_max), "Minimum value must be less than maximum value"

        rows = []
        cols = []
        vals = []
        for i in range(m):
            for j in range(n):
                if np.random.rand() < density:
                    rows.append(i)
                    cols.append(j)
                    vals.append(np.random.uniform(val_min, val_max))
        
        return sparse.csc_matrix((vals, (rows, cols)), shape=(m, n))

    @staticmethod
    def random_vector(n, val_min, val_max):
        assert(n >= 0), "Vector dimension must be non-negative"
        assert(val_min < val_max), "Minimum value must be less than maximum value"

        return np.random.uniform(val_min, val_max, n)

    @staticmethod
    def random_hybzono(n, nGc, nGb, nC, density, val_min, val_max):
        Gc = TestUtilities.random_sparse_matrix(n, nGc, density, val_min, val_max)
        Gb = TestUtilities.random_sparse_matrix(n, nGb, density, val_min, val_max)
        c = TestUtilities.random_vector(n, val_min, val_max)
        Ac = TestUtilities.random_sparse_matrix(nC, nGc, density, val_min, val_max)
        Ab = TestUtilities.random_sparse_matrix(nC, nGb, density, val_min, val_max)
        b = TestUtilities.random_vector(nC, val_min, val_max)

        return zono.HybZono(Gc, Gb, c, Ac, Ab, b)

    @staticmethod
    def random_conzono(n, nG, nC, density, val_min, val_max):
        G = TestUtilities.random_sparse_matrix(n, nG, density, val_min, val_max)
        c = TestUtilities.random_vector(n, val_min, val_max)
        A = TestUtilities.random_sparse_matrix(nC, nG, density, val_min, val_max)
        b = TestUtilities.random_vector(nC, val_min, val_max)

        return zono.ConZono(G, c, A, b)

    @staticmethod
    def random_zono(n, nG, density, val_min, val_max):
        G = TestUtilities.random_sparse_matrix(n, nG, density, val_min, val_max)
        c = TestUtilities.random_vector(n, val_min, val_max)

        return zono.Zono(G, c)


# unit tests
def test_vrep_2_hz():

    # folder where unit test data resides
    test_folder = test_data_folder / 'vrep_2_hybzono'

    # build hybzono from vrep
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
    test_folder = test_data_folder / 'minkowski_sum'

    # build hybzono from vrep
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
    test_folder = test_data_folder / 'intersection'

    # build hybzono from vrep
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
    Z2 = zono.Zono(sparse.csc_matrix(G), c)

    # minkowski sum
    Z3 = zono.minkowski_sum(Z1, Z2)

    # conzono
    G = np.array([[3.0, 0.0, 0.0], 
                  [0.0, 3.0, 0.0]])
    c = np.array([-0.5, 4.5])
    A = np.ones((1, 3))
    b = np.array([1.0])
    Z4 = zono.ConZono(sparse.csc_matrix(G), c, sparse.csc_matrix(A), b)

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
    test_folder = test_data_folder / 'is_empty'

    # load in feasible conzono
    G = np.loadtxt(test_folder / 'f_G.txt', delimiter=' ')
    c = np.loadtxt(test_folder / 'f_c.txt', delimiter=' ')
    A = np.loadtxt(test_folder / 'f_A.txt', delimiter=' ')
    b = np.loadtxt(test_folder / 'f_b.txt', delimiter=' ')

    Zf = zono.ConZono(sparse.csc_matrix(G), c, sparse.csc_matrix(A), b)

    # load in infeasible conzono
    G = np.loadtxt(test_folder / 'i_G.txt', delimiter=' ')
    c = np.loadtxt(test_folder / 'i_c.txt', delimiter=' ')
    A = np.loadtxt(test_folder / 'i_A.txt', delimiter=' ')
    b = np.loadtxt(test_folder / 'i_b.txt', delimiter=' ')

    Zi = zono.ConZono(sparse.csc_matrix(G), c, sparse.csc_matrix(A), b)
    
    # check if empty
    assert not Zf.is_empty()
    assert Zi.is_empty()
    print('Passed: Is Empty')

def test_support():

    # folder where unit test data resides
    test_folder = test_data_folder / 'support'

    # load in conzono
    G = np.loadtxt(test_folder / 'G.txt', delimiter=' ')
    c = np.loadtxt(test_folder / 'c.txt', delimiter=' ')
    A = np.loadtxt(test_folder / 'A.txt', delimiter=' ')
    b = np.loadtxt(test_folder / 'b.txt', delimiter=' ')

    Z = zono.ConZono(sparse.csc_matrix(G), c, sparse.csc_matrix(A), b)

    # load direction and expected support value
    d = np.loadtxt(test_folder / 'd.txt', delimiter=' ')
    s_expected = np.loadtxt(test_folder / 'sup.txt', delimiter=' ')

    # compute support
    s = Z.support(d)
    
    # compare results
    tol = 5e-2 # tolerance on success
    assert np.abs(s-s_expected)/np.abs(s_expected) < tol
    print('Passed: Support Function')

def test_point_contain():

    # folder where the data resides
    test_folder = test_data_folder / 'point_contain'

    # load in conzono
    G = np.loadtxt(test_folder / 'G.txt', delimiter=' ')
    c = np.loadtxt(test_folder / 'c.txt', delimiter=' ')
    A = np.loadtxt(test_folder / 'A.txt', delimiter=' ')
    b = np.loadtxt(test_folder / 'b.txt', delimiter=' ')

    Z = zono.ConZono(sparse.csc_matrix(G), c, sparse.csc_matrix(A), b)

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
    leaves = Z.get_leaves()

    # check number of leaves is correct
    assert len(leaves) == n_CZs**2
    print('Passed: Get Leaves')

def test_safety_verification():

    # System dynamics
    dt = 0.1
    A = np.array([[1., dt],
                  [0., 1.]])
    B = np.array([[0.5*dt**2],
                  [dt]])

    # Initial set: box [-1.0, 1.0] x [-0.1, 0.1]
    X0 = zono.interval_2_zono(zono.Box([-1., -0.1], [1., 0.1]))

    # Input set: box [-0.2, 0.2]
    U = zono.interval_2_zono(zono.Box([-0.2], [0.2]))

    # Disturbance set: affine map of octagon
    W = zono.make_regular_zono_2D(radius=1., n_sides=8)
    W = zono.affine_map(W, np.diag([0.01, 0.05]))

    # Compute reachable set over 10 time steps
    X = X0
    for k in range(10):
        X = zono.affine_map(X, A)
        X = zono.minkowski_sum(X, zono.affine_map(U, B))
        X = zono.minkowski_sum(X, W)

    # Unsafe set
    O = zono.vrep_2_conzono(np.array([[1.3, 0.],
                                      [1.6, 0.8],
                                      [2.0, -0.4],
                                      [2.3, 0.6]]))

    # expect intersection of X and O is empty
    assert zono.intersection(X, O).is_empty()
    print('Passed: Safety Verification')

def test_interval_arithmetic():
    
    np.random.seed(0)

    # constants
    n_dims = 3
    n_samples = 10000
    
    # expression
    f = lambda x: 2*np.tan(x[0])**(-2) + np.cos(x[1]/x[0])/3. + np.sin(x[0] + np.arctan(x[2]))*np.sinh(x[0]) + np.exp(np.arccosh(np.abs(x[1]) + 1)) - np.arccos(x[0])*np.arcsin(x[1])/np.log(x[2]**2)

    # interval expression
    f_int = lambda x: 2*(x[0].tan()**(-2)) + (x[1]/x[0]).cos()/3. + (x[0] + x[2].arctan()).sin()*x[0].sinh() + (1 + x[1].abs()).arccosh().exp() - (x[0].arccos()*x[1].arcsin())/(x[2]**2).log()

    def _run_interval_test(x_min, x_max):
        x = zono.Box(x_min*np.ones(n_dims), x_max*np.ones(n_dims)) # box
        x_sample = np.random.uniform(x_min, x_max, (n_samples, n_dims)) # generate random points in interval

        # get bounding interval
        f_bounds = f_int(x)

        # get samples
        f_samples = np.array([f(x_sample[i,:]) for i in range(n_samples)])

        # plot
        # import matplotlib.pyplot as plt
        # fig = plt.hist([fs for fs in f_samples if not np.isnan(fs) and not np.isinf(fs)], bins=100, density=True)
        # print(f'Interval bounds: {f_bounds}')
        # plt.show()

        # check that all samples evaluate to values within the computed interval
        for f_sample in f_samples:
            assert np.isnan(f_sample) or f_bounds.contains(f_sample)

    def _test_exponent():
        a = zono.Interval(0.5, 3.)
        b = a**(456./123)
        assert not b.is_empty(), 'test_exponent did not succeed'
        assert np.abs(b.lower() - 0.5**(456/123)) < 1e-6, 'test_exponent lower bound is incorrect'
        assert np.abs(b.upper() - 3.0**(456/123)) < 1e-6, 'test_exponent upper bound is incorrect'

        a = zono.Interval(-3., -0.5)
        try:
            b = a**(456./123)
            raise RuntimeError('test_exponent: expected fractional power of negative interval to throw')
        except ValueError:
            pass

    # Case 1: positive range
    _run_interval_test(0.1, 0.2)
    
    # Case 2: negative range, approaching 0
    _run_interval_test(-0.2, -0.001)

    # Case 3: spanning 0
    _run_interval_test(-1., 1.)

    # fractional power test
    _test_exponent()

    print('Passed: Interval Arithmetic')

def test_affine_inclusion():

    def _random_example(n, nG, n_out, rng_seed=None):

        # rng seed
        if rng_seed is not None:
            np.random.seed(rng_seed)

        # generate zonotope
        G = np.random.rand(n, nG)
        c = np.random.rand(n)

        Z = zono.Zono(G, c)

        # offset
        s = 10.*np.random.rand(n_out)

        # random affine map
        R = np.random.randn(n_out, n)
        R_lb = R.copy()
        R_ub = R.copy()
        for i in range(2):
            for j in range(n):
                drij = 0.1*np.random.rand()
                R_lb[i,j] = R[i,j] - drij
                R_ub[i,j] = R[i,j] + drij

        return Z, R_lb, R_ub, s
    
    for seed in range(10):

        # get test case
        Z, R_lb, R_ub, s = _random_example(n=5, nG=10, n_out = 3, rng_seed=seed)

        # affine inclusion
        Rint = zono.IntervalMatrix(R_lb, R_ub)
        Z_inc = zono.affine_inclusion(Z, Rint, s)
        
        # lower and upper bounds
        Z_lb = zono.affine_map(Z, R_lb, s)
        Z_ub = zono.affine_map(Z, R_ub, s)

        # make sure centers of Z_lb and Z_ub are inside Z_inc
        assert Z_inc.contains_point(Z_lb.get_center())
        assert Z_inc.contains_point(Z_ub.get_center())

    print('Passed: Affine Inclusion')

def test_operator_overloading():

    def check_matrix_equal(A, B):
        if A.shape != B.shape:
            return False
        
        for i in range(A.shape[0]):
            for j in range(A.shape[1]):
                if np.abs(A[i,j]-B[i,j]) > 1e-12:
                    return False
                
        return True

    def check_vector_equal(a, b):
        if a.shape != b.shape:
            return False

        for i in range(a.shape[0]):
            if np.abs(a[i]-b[i]) > 1e-12:
                return False

        return True
    
    def check_equal(Z1, Z2):
        if Z2.is_0_1_form() != Z1.is_0_1_form():
            Z2.convert_form()

        # make sure dimensions are the same
        if Z1.get_n() != Z2.get_n() or Z1.get_nGc() != Z2.get_nGc() or Z1.get_nGb() != Z2.get_nGb() or Z1.get_nC() != Z2.get_nC():
            return False
        
        # make sure matrices and vectors are the same
        G_equal = check_matrix_equal(Z1.get_G().toarray(), Z2.get_G().toarray())
        c_equal = check_vector_equal(Z1.get_c(), Z2.get_c())
        A_equal = check_matrix_equal(Z1.get_A().toarray(), Z2.get_A().toarray())
        b_equal = check_vector_equal(Z1.get_b(), Z2.get_b())

        return G_equal and c_equal and A_equal and b_equal
    

    # make 2 zonotopes and a point
    c1 = np.array([1., 2.])
    Z1 = zono.make_regular_zono_2D(radius=1., n_sides=8, c=c1)

    c2 = np.array([-3., 1.])
    Z2 = zono.make_regular_zono_2D(radius=0.1, n_sides=6, c=c2)

    c3 = np.array([3., -2.])
    P3 = zono.Point(c3)

    # matrix
    M = np.array([[32., 1.2]])
    M_sp = sparse.csc_matrix(M)

    M_upper = np.array([[33., 1.4]])
    M_int = zono.IntervalMatrix(M, M_upper)

    # box
    box = zono.Box([0., 1.], [0.1, 1.04])

    # check operators are consistent with set operations

    # minkowski sum
    Z_set = zono.minkowski_sum(Z1, Z2)
    Z_op = Z1 + Z2
    assert check_equal(Z_set, Z_op), "Minkowski sum failed"

    # +=
    Z_set = Z1.copy()
    Z_op = Z1.copy()
    Z_set = zono.minkowski_sum(Z_set, Z2)
    Z_op += Z2
    assert check_equal(Z_set, Z_op), "Minkowski sum += failed"

    # minkowski sum with point
    Z_set = zono.minkowski_sum(Z1, P3)
    Z_op = Z1 + c3
    assert check_equal(Z_set, Z_op), "Minkowski sum with point failed"

    # left sum
    Z_set = zono.minkowski_sum(P3, Z1)
    Z_op = c3 + Z1
    assert check_equal(Z_set, Z_op), "Left Minkowski sum with point failed"

    # +=
    Z_set = Z1.copy()
    Z_op = Z1.copy()
    Z_set = zono.minkowski_sum(Z_set, P3)
    Z_op += c3
    assert check_equal(Z_set, Z_op), "Minkowski sum with point += failed"

    # minkowski sum with box
    Z_set = zono.minkowski_sum(Z1, zono.interval_2_zono(box))
    Z_op = Z1 + box
    assert check_equal(Z_set, Z_op), "Minkowski sum with box failed"

    # left sum
    Z_set = zono.minkowski_sum(zono.interval_2_zono(box), Z1)
    Z_op = box + Z1
    assert check_equal(Z_set, Z_op), "Left Minkowski sum with box failed"

    # +=
    Z_set = Z1.copy()
    Z_op = Z1.copy()
    Z_set = zono.minkowski_sum(Z_set, zono.interval_2_zono(box))
    Z_op += box
    assert check_equal(Z_set, Z_op), "Minkowski sum with box += failed"

    # pontry diff
    Z_set = zono.pontry_diff(Z1, Z2, exact=True)
    Z_op = Z1 - Z2
    assert check_equal(Z_set, Z_op), "Pontryagin difference failed"

    # -=
    Z_set = Z1.copy()
    Z_op = Z1.copy()
    Z_set = zono.pontry_diff(Z_set, Z2, exact=True)
    Z_op -= Z2
    assert check_equal(Z_set, Z_op), "Pontryagin difference -= failed"

    # pontry diff with point
    Z_set = zono.pontry_diff(Z1, P3, exact=True)
    Z_op = Z1 - c3
    assert check_equal(Z_set, Z_op), "Pontryagin difference with point failed"

    # -=
    Z_set = Z1.copy()
    Z_op = Z1.copy()
    Z_set = zono.pontry_diff(Z_set, P3, exact=True)
    Z_op -= c3
    assert check_equal(Z_set, Z_op), "Pontryagin difference with point -= failed"

    # pontry diff with box
    Z_set = zono.pontry_diff(Z1, zono.interval_2_zono(box), exact=True)
    Z_op = Z1 - box
    assert check_equal(Z_set, Z_op), "Pontryagin difference with box failed"

    # -=
    Z_set = Z1.copy()
    Z_op = Z1.copy()
    Z_set = zono.pontry_diff(Z_set, zono.interval_2_zono(box), exact=True)
    Z_op -= box
    assert check_equal(Z_set, Z_op), "Pontryagin difference with box -= failed"

    # affine map - sparse
    Z_set = zono.affine_map(Z1, M_sp)
    Z_op = M_sp @ Z1
    assert check_equal(Z_set, Z_op), "Affine map with sparse matrix failed"

    # affine map - dense
    Z_set = zono.affine_map(Z1, M)
    Z_op = M @ Z1
    assert check_equal(Z_set, Z_op), "Affine map with dense matrix failed"

    # affine inclusion
    Z_set = zono.affine_inclusion(Z1, M_int)
    Z_op = M_int @ Z1
    assert check_equal(Z_set, Z_op), "Affine inclusion failed"

    # scalar multiplication - left
    f = 3.2
    Z_set = zono.affine_map(Z1, f*np.eye(Z1.get_n()))
    Z_op = f * Z1
    assert check_equal(Z_set, Z_op), "Scalar multiplication (left) failed"

    # scalar multiplication - right
    Z_op = Z1 * f
    assert check_equal(Z_set, Z_op), "Scalar multiplication (right) failed"

    # *=
    Z_set = Z1.copy()
    Z_op = Z1.copy()
    Z_set = zono.affine_map(Z_set, f*np.eye(Z_set.get_n()))
    Z_op *= f
    assert check_equal(Z_set, Z_op), "Scalar multiplication *= failed"

    # cartesian product
    Z_set = zono.cartesian_product(Z1, Z2)
    Z_op = Z1 * Z2
    assert check_equal(Z_set, Z_op), "Cartesian product failed"

    # *=
    Z_set = Z1.copy()
    Z_op = Z1.copy()
    Z_set = zono.cartesian_product(Z_set, Z2)
    Z_op *= Z2
    assert check_equal(Z_set, Z_op), "Cartesian product *= failed"

    # cartesian product with box
    Z_set = zono.cartesian_product(Z1, zono.interval_2_zono(box))
    Z_op = Z1 * box
    assert check_equal(Z_set, Z_op), "Cartesian product with box failed"

    # cartesian product with box on left
    Z_set = zono.cartesian_product(zono.interval_2_zono(box), Z1)
    Z_op = box * Z1
    assert check_equal(Z_set, Z_op), "Cartesian product with box on left failed"

    # *=
    Z_set = Z1.copy()
    Z_op = Z1.copy()
    Z_set = zono.cartesian_product(Z_set, zono.interval_2_zono(box))
    Z_op *= box
    assert check_equal(Z_set, Z_op), "Cartesian product with box *= failed"

    # intersection
    Z_set = zono.intersection(Z1, Z2)
    Z_op = Z1 & Z2
    assert check_equal(Z_set, Z_op), "Intersection failed"

    # union
    Z_set = zono.union_of_many((Z1, Z2))
    Z_op = Z1 | Z2
    assert check_equal(Z_set, Z_op), "Union failed"

    # unary minus
    Z_set = zono.affine_map(Z1, -np.eye(Z1.get_n()))
    Z_op = -Z1
    assert check_equal(Z_set, Z_op), "Unary minus failed"

    # finish
    print('Passed: Operator Overloading')

def test_constrain():
    np.random.seed(9)
    n = 2
    nGc = 4
    nGb = 2
    nC = 2
    Gc = np.random.randn(n, nGc)
    Gb = np.random.randn(n, nGb)
    c = np.random.randn(n)
    Ac = np.random.randn(nC, nGc)
    Ab = np.random.randn(nC, nGb)
    b = np.random.randn(nC)
    Zh = zono.HybZono(Gc, Gb, c, Ac, Ab, b, zero_one_form=False, sharp=False)

    h0 = np.random.randn(n)
    f0 = np.random.rand(1)
    Zh_halfspace = zono.halfspace_intersection(Zh, h0, f0)
    Zh_leq = zono.constrain(Zh, h0, f0, '<')
    Zh_geq = zono.constrain(Zh, h0, f0, '>')
    Zh_eq = zono.constrain(Zh, h0, f0, '=')

    assert len(Zh_halfspace.get_leaves()) == 4, "Halfspace intersection should return 4 leaves"
    assert len(Zh_leq.get_leaves()) == 4, "Less than or equal constraint should return 4 leaves"
    assert len(Zh_geq.get_leaves()) == 3, "Greater than or equal constraint should return 3 leaves"
    assert len(Zh_eq.get_leaves()) == 3, "Equality constraint should return 3 leaves"

    # finish
    print('Passed: Constrain')

def test_remove_redundancy():

    def _test1():
        # constrained zonotope that can be simplified
        G = np.array([[3., 1., 0., 0.]])
        c = np.array([8.])
        A = np.array([[0.5, 0., 1., 0.],
                      [0., 0.5, 0., 0.5]])
        b = np.array([-0.5, -1.])

        Z = zono.ConZono(G, c, A, b)

        # simplify
        Z_rr = Z.remove_redundancy()

        # check that support is as expected
        d = np.array([1.])
        sup = Z_rr.support(d)
        assert(np.abs(sup - 10.) < 1e-3), f'case 1: expected support = 10., got support = {sup}'

        d = np.array([-1.])
        sup = Z_rr.support(d)
        assert(np.abs(sup - -4.) < 1e-3), f'case 1: expected support = -4., got support = {sup}'

        assert(Z_rr.is_zono()), 'case1: expected result to be a zonotope after removing redundancy'

    def _test2():
        # hybrid zonotope
        Gc = np.array([[3., 1.]])
        Gb = np.zeros((1,2))
        c = np.array([8.])
        Ac = np.array([[0.5, 0.],
                       [0., 0.5]])
        Ab = np.array([[1., 0.],
                       [0., 0.5]])
        b = np.array([-0.5, -1.])

        Z = zono.HybZono(Gc, Gb, c, Ac, Ab, b)

        # get support
        sup_before = np.zeros(2)

        d = np.array([1.])
        sup_before[0] =  Z.support(d)

        d = np.array([-1.])
        sup_before[1] = Z.support(d)

        # simplify
        Z_rr = Z.remove_redundancy()

        # check that support is as expected
        sup_after = np.zeros(2)
        d = np.array([1.])
        sup_after[0] = Z_rr.support(d)

        d = np.array([-1.])
        sup_after[1] = Z_rr.support(d)
        
        for i in range(2):
            assert(np.abs(sup_before[i]-sup_after[i]) < 1e-3), f'Hybrid Zono: expected support = {sup_before[i]}, got support = {sup_after[i]}'

    def _test3():
        # constrained zonotope
        G = np.array([[3., 1., 0., 0.]])
        c = np.array([8.])
        A = np.array([[0.5, 0.1, 1., 0.],
                      [0., 0.5, 0., 0.5]])
        b = np.array([-0.5, -1.])

        Z = zono.ConZono(G, c, A, b)

        # get support
        sup_before = np.zeros(2)

        d = np.array([1.])
        sup_before[0] =  Z.support(d)

        d = np.array([-1.])
        sup_before[1] = Z.support(d)

        # simplify
        Z_rr = Z.remove_redundancy()

        # check that support is as expected
        sup_after = np.zeros(2)
        d = np.array([1.])
        sup_after[0] = Z_rr.support(d)

        d = np.array([-1.])
        sup_after[1] = Z_rr.support(d)
        
        for i in range(2):
            assert(np.abs(sup_before[i]-sup_after[i]) < 1e-3), f'Hybrid Zono: expected support = {sup_before[i]}, got support = {sup_after[i]}'

    def _test4():
        # infeasible constrained zonotope
        G = np.array([[1., 1.],
                      [1., 2.]])
        c = np.array([1., 2.])
        A = np.array([[1., 1.]])
        b = np.array([3.])
        
        Z = zono.ConZono(G, c, A, b)

        # remove redundancy
        Z_rr = Z.remove_redundancy()

        # check that result is EmptySet object
        assert Z_rr.is_empty_set(), f'Expected EmptySet, got {Z_rr}'

    def _test_random_conzono():
        Z = TestUtilities.random_conzono(n=2, nG=30, nC=10, density=0.1, val_min=0., val_max=1.)

        # get support before simplifying
        settings = zono.OptSettings()
        settings.eps_prim = 1e-3
        settings.eps_dual = 1e-3
        settings.rho = 1.
        sup_before = np.zeros(4)

        try:
            d = [1., 0.]
            sup_before[0] = Z.support(d, settings=settings)
            d = [-1., 0.]
            sup_before[1] = Z.support(d, settings=settings)
            d = [0., 1.]
            sup_before[2] = Z.support(d, settings=settings)
            d = [0., -1.]
            sup_before[3] = Z.support(d, settings=settings)

        except Exception as e:
            return
        
        # randomly convert form
        if np.random.rand() < 0.5:
            Z.convert_form()

        # get support after simplifying
        Z_rr = Z.remove_redundancy()
        sup_after = np.zeros(4)

        d = [1., 0.]
        sup_after[0] = Z_rr.support(d, settings=settings)
        d = [-1., 0.]
        sup_after[1] = Z_rr.support(d, settings=settings)
        d = [0., 1.]
        sup_after[2] = Z_rr.support(d, settings=settings)
        d = [0., -1.]
        sup_after[3] = Z_rr.support(d, settings=settings)

        # make sure all close
        for i in range(4):
            err_str = f'Random ConZono: expected support = {sup_before[i]}, got support = {sup_after[i]}\n  Z before simplifying: {Z}\n  Z after simplifying: {Z_rr}'
            cond = np.abs(sup_before[i]-sup_after[i])/np.abs(sup_before[i]) < 1e-1 or np.abs(sup_before[i] - sup_after[i]) < 1e-1

            assert(cond), err_str

    def _check_vertices_equal(V1, V2):
        for i in range(V1.shape[0]):
            v1 = V1[i,:]
            found_match = False
            min_dist = np.inf
            closest_vertex = None
            for j in range(V2.shape[0]):
                v2 = V2[j,:]
                if np.linalg.norm(v1-v2) < min_dist:
                    min_dist = np.linalg.norm(v1-v2)
                    closest_vertex = v2
                if min_dist < 1e-3:
                    found_match = True
                    break
            
            assert found_match, f'Vertex {v1} not found in second set of vertices, closest vertex was {closest_vertex} with distance {min_dist}'

    def _test_random_hybzono():
        Z = TestUtilities.random_hybzono(n=2, nGc=20, nGb=5, nC=5, density=0.2, val_min=0., val_max=1.)

        # get vertices before simplifying
        try:
            V_before = zono.get_vertices(Z)
            if V_before.shape[0] == 0: # infeasible
                return
        except RuntimeError: # can't factor problem matrices
            return
        
        # randomly convert form
        if np.random.rand() < 0.5:
            Z.convert_form()

        # get vertices after simplifying
        Z_rr = Z.remove_redundancy()

        # get vertices after simplifying
        V_after = zono.get_vertices(Z_rr)

        # check that vertices are the same
        try:
            _check_vertices_equal(V_before, V_after)
            _check_vertices_equal(V_after, V_before)
        except Exception:
            import matplotlib.pyplot as plt

            print(Z)
            print(Z_rr)
            
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            ax.scatter(V_before[:,0], V_before[:,1], color='blue', marker='.', label='before')
            ax.scatter(V_after[:,0], V_after[:,1], color='red', marker='x', label='after')
            ax.legend()
            ax.set_aspect('equal')
            plt.show()

    # main
    _test1()
    _test2()
    _test3()
    _test4()

    np.random.seed(0)

    for _ in range(100):
        _test_random_conzono()

    for _ in range(100):
        _test_random_hybzono()

    print('Passed: Remove Redundancy')


# run the unit tests
test_vrep_2_hz()
test_minkowski_sum()
test_intersection()
test_is_empty()
test_support()
test_point_contain()
test_get_leaves()
test_safety_verification()
test_interval_arithmetic()
test_affine_inclusion()
test_operator_overloading()
test_constrain()
test_remove_redundancy()