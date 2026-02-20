import collections.abc
import flags
import numpy
import numpy.typing
import scipy.sparse
import typing
from typing import Any, ClassVar, overload

EQUAL: IneqType
GREATER: IneqType
GREATER_OR_EQUAL: IneqType
LESS: IneqType
LESS_OR_EQUAL: IneqType

class Box:
    """Box (i.e., interval vector) class"""
    def __init__(self, x_lb: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], x_ub: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, '[m, 1]']) -> None:
        '''__init__(self: zonoopt._core.Box, x_lb: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], x_ub: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"]) -> None


                        Constructor from intervals of lower and upper bounds

                        Args:
                            x_lb (numpy.array): vector of lower bounds
                            x_ub (numpy.array): vector of upper bounds
            
        '''
    def center(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], '[m, 1]']:
        '''center(self: zonoopt._core.Box) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]


                        Gets center of box (x_ub + x_lb) / 2

                        Returns:
                            numpy.array: center of interval
            
        '''
    def contract(self, A: scipy.sparse.csr_matrix[numpy.float64], b: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], iter: typing.SupportsInt | typing.SupportsIndex) -> bool:
        '''contract(self: zonoopt._core.Box, A: scipy.sparse.csr_matrix[numpy.float64], b: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], iter: typing.SupportsInt | typing.SupportsIndex) -> bool


                        Interval contractor.

                        Executes a forward-backward interval contractor for the equality constraint A*x=b.
                        For points x in the box, this shrinks the box without removing any points x that satisfy A*x=b.
                        If the contractor detects that the box does not intersect A*x=b, then this function will return false.

                        Args:
                            A (scipy.sparse.csr_matrix): constraint matrix
                            b (numpy.vector): constraint vector
                            iter (int): number of contractor iterations

                        Returns:
                            bool: flag indicating that the contractor did not detect that A*x=b and the box do not intersect
            
        '''
    def copy(self) -> Box:
        """copy(self: zonoopt._core.Box) -> zonoopt._core.Box


                        Copies Box object

                        Returns:
                            Box: copy of object
            
        """
    def dot(self, x: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, '[m, 1]']) -> Interval:
        '''dot(self: zonoopt._core.Box, x: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"]) -> zonoopt._core.Interval


                        Linear map with vector

                        Args:
                            x (numpy.array): vector

                        Returns:
                            Interval: result of linear map of box with vector
            
        '''
    def linear_map(self, A: scipy.sparse.csr_matrix[numpy.float64]) -> Box:
        """linear_map(self: zonoopt._core.Box, A: scipy.sparse.csr_matrix[numpy.float64]) -> zonoopt._core.Box


                        Linear map of box based on interval arithmetic

                        Args:
                            A (scipy.sparse.csr_matrix): linear map matrix

                        Returns:
                            Box: linear mapped box
            
        """
    def lower(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], '[m, 1]']:
        '''lower(self: zonoopt._core.Box) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]


                        Get reference to lower bounds

                        Returns:
                            numpy.array: lower bounds
            
        '''
    def project(self, x: typing.Annotated[numpy.typing.NDArray[numpy.float64], '[m, 1]', 'flags.writeable']) -> None:
        '''project(self: zonoopt._core.Box, x: typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]", "flags.writeable"]) -> None


                        Projects vector onto the Box (in place)

                        Args:
                            x (numpy.array): vector to be projected
            
        '''
    def radius(self) -> Box:
        """radius(self: zonoopt._core.Box) -> zonoopt._core.Box


                        Get radius of box

                        Returns box with intervals centered at zero with width equal to the width of the original box

                        Returns:
                            Box: radius of box
            
        """
    def size(self) -> int:
        """size(self: zonoopt._core.Box) -> int


                        Get size of Box object

                        Returns:
                            int: size of box
            
        """
    def upper(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], '[m, 1]']:
        '''upper(self: zonoopt._core.Box) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]


                        Get reference to upper bounds

                        Returns:
                            numpy.array: upper bounds
            
        '''
    def width(self) -> float:
        """width(self: zonoopt._core.Box) -> float


                        Get width of box.

                        Specifically, this returns the max width for any interval in the box

                        Returns:
                            float: width of box
            
        """
    def __add__(self, other: Box) -> Box:
        """__add__(self: zonoopt._core.Box, other: zonoopt._core.Box) -> zonoopt._core.Box


                        Elementwise addition

                        Args:
                            other (Box): rhs box

                        Returns:
                            Box: self + other (elementwise)
            
        """
    def __getitem__(self, i: typing.SupportsInt | typing.SupportsIndex) -> Interval:
        """__getitem__(self: zonoopt._core.Box, i: typing.SupportsInt | typing.SupportsIndex) -> zonoopt._core.Interval


                        Get interval at index i

                        Args:
                            i (int): index

                        Returns:
                            Interval: interval at index i in Box
            
        """
    @overload
    def __mul__(self, other: Box) -> Box:
        """__mul__(*args, **kwargs)
        Overloaded function.

        1. __mul__(self: zonoopt._core.Box, other: zonoopt._core.Box) -> zonoopt._core.Box


                        Elementwise multiplication

                        Args:
                            other (Box): rhs box

                        Returns:
                            Box: self * other (elementwise)
            

        2. __mul__(self: zonoopt._core.Box, alpha: typing.SupportsFloat | typing.SupportsIndex) -> zonoopt._core.Box


                        Elementwise multiplication with scalar

                        Args:
                            alpha (float): scalar multiplier

                        Returns:
                            Box: alpha * self (elementwise)
            
        """
    @overload
    def __mul__(self, alpha: typing.SupportsFloat | typing.SupportsIndex) -> Box:
        """__mul__(*args, **kwargs)
        Overloaded function.

        1. __mul__(self: zonoopt._core.Box, other: zonoopt._core.Box) -> zonoopt._core.Box


                        Elementwise multiplication

                        Args:
                            other (Box): rhs box

                        Returns:
                            Box: self * other (elementwise)
            

        2. __mul__(self: zonoopt._core.Box, alpha: typing.SupportsFloat | typing.SupportsIndex) -> zonoopt._core.Box


                        Elementwise multiplication with scalar

                        Args:
                            alpha (float): scalar multiplier

                        Returns:
                            Box: alpha * self (elementwise)
            
        """
    def __setitem__(self, i: typing.SupportsInt | typing.SupportsIndex, val: Interval) -> None:
        """__setitem__(self: zonoopt._core.Box, i: typing.SupportsInt | typing.SupportsIndex, val: zonoopt._core.Interval) -> None


                        Set indexed interval in box to specified value

                        Args:
                            i (int): index
                            val (Interval): new interval for index i in Box
            
        """
    def __sub__(self, other: Box) -> Box:
        """__sub__(self: zonoopt._core.Box, other: zonoopt._core.Box) -> zonoopt._core.Box


                        Elementwise subtraction

                        Args:
                            other (Box): rhs box

                        Returns:
                            Box: self - other (elementwise)
            
        """
    def __truediv__(self, other: Box) -> Box:
        """__truediv__(self: zonoopt._core.Box, other: zonoopt._core.Box) -> zonoopt._core.Box


                        Elementwise division

                        Args:
                            other (Box): rhs box

                        Returns:
                            Box: self / other (elementwise)
            
        """

class ConZono(HybZono):
    """
                Constrained zonotope class
                
                A constrained zonotope is defined as:
                Z = {G \\xi + c | A \\xi = b, \\xi in [-1, 1]^nG}.
                Equivalently, the following shorthand can be used: Z = <G, c, A, b>.
                Optionally, in 0-1 form, the factors are xi in [0,1].
                The set dimension is n, and the number of equality constraints is nC.
            """
    def __init__(self, G: scipy.sparse.csc_matrix[numpy.float64], c: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], A: scipy.sparse.csc_matrix[numpy.float64], b: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], zero_one_form: bool = ...) -> None:
        '''__init__(self: zonoopt._core.ConZono, G: scipy.sparse.csc_matrix[numpy.float64], c: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], A: scipy.sparse.csc_matrix[numpy.float64], b: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], zero_one_form: bool = False) -> None


                        ConZono constructor
                
                        Args:
                            G (scipy.sparse.csc_matrix): generator matrix
                            c (numpy.array): center
                            A (scipy.sparse.csc_matrix): constraint matrix
                            b (numpy.array): constraint vector
                            zero_one_form (bool, optional): true if set is in 0-1 form
            
        '''
    def constraint_reduction(self) -> None:
        """constraint_reduction(self: zonoopt._core.ConZono) -> None


                        Execute constraint reduction algorithm from Scott et. al. 2016

                        Removes one constraint and one generator from the constrained zonotope.
                        The resulting set is an over-approximation of the original set.
            
        """
    def set(self, G: scipy.sparse.csc_matrix[numpy.float64], c: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], A: scipy.sparse.csc_matrix[numpy.float64], b: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], zero_one_form: bool = ...) -> None:
        '''set(self: zonoopt._core.ConZono, G: scipy.sparse.csc_matrix[numpy.float64], c: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], A: scipy.sparse.csc_matrix[numpy.float64], b: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], zero_one_form: bool = False) -> None


                        Reset constrained zonotope object with the given parameters.
                
                        Args:
                            G (scipy.sparse.csc_matrix): generator matrix
                            c (numpy.array): center
                            A (scipy.sparse.csc_matrix): constraint matrix
                            b (numpy.array): constraint vector
                            zero_one_form (bool, optional): true if set is in 0-1 form
            
        '''
    def to_zono_approx(self, *args, **kwargs):
        """to_zono_approx(self: zonoopt._core.ConZono) -> ZonoOpt::Zono


                        Compute outer approximation of constrained zonotope as zonotope using SVD

                        Returns:
                            Zono: Zonotope over-approximation
            
        """

class EmptySet(ConZono):
    """
                Empty Set class

                Used to facilitate set operations with trivial solutions when one of the sets is an empty set.
            """
    def __init__(self, n: typing.SupportsInt | typing.SupportsIndex) -> None:
        """__init__(self: zonoopt._core.EmptySet, n: typing.SupportsInt | typing.SupportsIndex) -> None


                        EmptySet constructor

                        Args:
                            n (int): dimension
            
        """

class HybZono:
    """
            Hybrid zonotope class
             
            A hybrid zonotope is defined as:
            Z = {Gc * xi_c + Gb * xi_b + c | Ac * xi_c + Ab * xi_b = b, xi_c in [-1, 1]^nGc, xi_b in {-1, 1}^nGb}.
            Equivalently, the following shorthand can be used: Z = <Gc, Gb, c, Ac, Ab, b>.
            Optionally, in 0-1 form, the factors are xi_c in [0, 1]^nGc, xi_b in {0, 1}^nGb. 
            The set dimension is n, and the number of equality constraints is nC.
        """
    def __init__(self, Gc: scipy.sparse.csc_matrix[numpy.float64], Gb: scipy.sparse.csc_matrix[numpy.float64], c: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], Ac: scipy.sparse.csc_matrix[numpy.float64], Ab: scipy.sparse.csc_matrix[numpy.float64], b: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], zero_one_form: bool = ..., sharp: bool = ...) -> None:
        '''__init__(self: zonoopt._core.HybZono, Gc: scipy.sparse.csc_matrix[numpy.float64], Gb: scipy.sparse.csc_matrix[numpy.float64], c: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], Ac: scipy.sparse.csc_matrix[numpy.float64], Ab: scipy.sparse.csc_matrix[numpy.float64], b: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], zero_one_form: bool = False, sharp: bool = False) -> None


                        HybZono constructor
                
                        Args:
                            Gc (scipy.sparse.csc_matrix): continuous generator matrix
                            Gb (scipy.sparse.csc_matrix): binary generator matrix
                            c (numpy.array): center
                            Ac (scipy.sparse.csc_matrix): continuous constraint matrix
                            Ab (scipy.sparse.csc_matrix): binary constraint matrix
                            b (numpy.array): constraint vector
                            zero_one_form (bool, optional): true if set is in 0-1 form
                            sharp (bool, optional): true if set is known to be sharp, i.e., convex relaxation = convex hull
            
        '''
    def bounding_box(self, settings: OptSettings = ..., solution: OptSolution = ..., warm_start_params: WarmStartParams = ...) -> Box:
        """bounding_box(self: zonoopt._core.HybZono, settings: zonoopt._core.OptSettings = OptSettings structure: verbose: false verbosity_interval: 100 t_max: 1.79769e+308 k_max_admm: 5000 rho: 10 eps_dual: 0.01 eps_prim: 0.001 k_inf_check: 10 inf_norm_conv: true use_interval_contractor: true contractor_iter: 1 search_mode: 0 polish: 1 eps_dual_search: 0.1 eps_prim_search: 0.01 eps_r: 0.01 eps_a: 0.1 k_max_bnb: 100000 n_threads_bnb: 4 n_threads_admm_fp: 3 single_threaded_admm_fp: false max_nodes: 100000 contractor_tree_search_depth: 10 enable_perturb_admm_fp: true k_max_admm_fp_ph1: 10000 k_max_admm_fp_ph2: 90000 cycle_detection_buffer_size: 20 eps_perturb: 0.001 k_restart: 5000 enable_rng_seed: false rng_seed: 0 enable_restart_admm_fp: true, solution: zonoopt._core.OptSolution = None, warm_start_params: zonoopt._core.WarmStartParams = <zonoopt._core.WarmStartParams object at 0x73c4f3897d70>) -> zonoopt._core.Box


                        Computes a bounding box of the set object as a Box object.
                
                        Args:
                            settings (OptSettings, optional): optimization settings structure
                            solution (OptSolution, optional): optimization solution structure pointer, populated with result
                            warm_start_params (WarmStartParams, optional): warm start parameters structure

                        Returns:
                            Box: bounding box of the set

                        In general, solves 2*n support optimizations where n is the set dimension to compute a bounding box.
            
        """
    def complement(self, delta_m: typing.SupportsFloat | typing.SupportsIndex = ..., remove_redundancy: bool = ..., settings: OptSettings = ..., solution: OptSolution = ..., n_leaves: typing.SupportsInt | typing.SupportsIndex = ..., contractor_iter: typing.SupportsInt | typing.SupportsIndex = ...) -> HybZono:
        '''complement(self: zonoopt._core.HybZono, delta_m: typing.SupportsFloat | typing.SupportsIndex = 100, remove_redundancy: bool = True, settings: zonoopt._core.OptSettings = OptSettings structure: verbose: false verbosity_interval: 100 t_max: 1.79769e+308 k_max_admm: 5000 rho: 10 eps_dual: 0.01 eps_prim: 0.001 k_inf_check: 10 inf_norm_conv: true use_interval_contractor: true contractor_iter: 1 search_mode: 0 polish: 1 eps_dual_search: 0.1 eps_prim_search: 0.01 eps_r: 0.01 eps_a: 0.1 k_max_bnb: 100000 n_threads_bnb: 4 n_threads_admm_fp: 3 single_threaded_admm_fp: false max_nodes: 100000 contractor_tree_search_depth: 10 enable_perturb_admm_fp: true k_max_admm_fp_ph1: 10000 k_max_admm_fp_ph2: 90000 cycle_detection_buffer_size: 20 eps_perturb: 0.001 k_restart: 5000 enable_rng_seed: false rng_seed: 0 enable_restart_admm_fp: true, solution: zonoopt._core.OptSolution = None, n_leaves: typing.SupportsInt | typing.SupportsIndex = 2147483647, contractor_iter: typing.SupportsInt | typing.SupportsIndex = 100) -> zonoopt._core.HybZono


                    Computes the complement of the set Z.
            
                    Args:
                        delta_m (float, optional): parameter defining range of complement
                        remove_redundancy (bool, optional): remove redundant constraints and unused generators in get_leaves function call
                        settings (OptSettings, optional): optimization settings for get_leaves function call
                        solution (OptSolution, optional): optimization solution for get_leaves function call
                        n_leaves (int, optional): maximum number of leaves to return in get_leaves function call
                        contractor_iter (int, optional): number of interval contractor iterations in remove_redundancy if using
            
                    Returns:
                        HybZono: Hybrid zonotope complement of the given set
            
                    Computes the complement according to the method of Bird and Jain:
                    "Unions and Complements of Hybrid Zonotopes"
                    delta_m is a parameter that defines the set over which the complement is defined.
                    For a constrained zonotope, the complement is restricted to the set
                    X = {G \\xi + c | A \\xi = b, \\xi \\in [-1-delta_m, 1+delta+m]^{nG}}.
            
        '''
    def contains_point(self, x: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], settings: OptSettings = ..., solution: OptSolution = ..., warm_start_params: WarmStartParams = ...) -> bool:
        '''contains_point(self: zonoopt._core.HybZono, x: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], settings: zonoopt._core.OptSettings = OptSettings structure: verbose: false verbosity_interval: 100 t_max: 1.79769e+308 k_max_admm: 5000 rho: 10 eps_dual: 0.01 eps_prim: 0.001 k_inf_check: 10 inf_norm_conv: true use_interval_contractor: true contractor_iter: 1 search_mode: 0 polish: 1 eps_dual_search: 0.1 eps_prim_search: 0.01 eps_r: 0.01 eps_a: 0.1 k_max_bnb: 100000 n_threads_bnb: 4 n_threads_admm_fp: 3 single_threaded_admm_fp: false max_nodes: 100000 contractor_tree_search_depth: 10 enable_perturb_admm_fp: true k_max_admm_fp_ph1: 10000 k_max_admm_fp_ph2: 90000 cycle_detection_buffer_size: 20 eps_perturb: 0.001 k_restart: 5000 enable_rng_seed: false rng_seed: 0 enable_restart_admm_fp: true, solution: zonoopt._core.OptSolution = None, warm_start_params: zonoopt._core.WarmStartParams = <zonoopt._core.WarmStartParams object at 0x73c4f38cd070>) -> bool


                        Checks whether the point x is contained in the set object.
                
                        Args:
                            x (numpy.array): point to be checked for set containment
                            settings (OptSettings, optional): optimization settings structure
                            solution (OptSolution, optional): optimization solution structure pointer, populated with result
                            warm_start_params (WarmStartParams, optional): warm start parameters structure

                        Returns:
                            bool: true if set contains point, false otherwise

                        False positives are possible; will return true if the optimization converges within the specified tolerances.
                        Will return false only if an infeasibility certificate is found, i.e., false negatives are not possible.
            
        '''
    def convert_form(self) -> None:
        """convert_form(self: zonoopt._core.HybZono) -> None


                        Converts the set representation between -1-1 and 0-1 forms.
                
                        This method converts the set representation between -1-1 and 0-1 forms. 
                        If the set is in -1-1 form, then xi_c in [-1,1] and xi_b in {-1,1}.
                        If the set is in 0-1 form, then xi_c in [0,1] and xi_b in {0,1}.
            
        """
    def convex_relaxation(self, *args, **kwargs):
        """convex_relaxation(self: zonoopt._core.HybZono) -> ZonoOpt::ConZono


                        Computes the convex relaxation of the hybrid zonotope.
                
                        Returns:
                            ConZono: Constrained zonotope Z = <[Gc, Gb], c, [Ac, Ab,], b>

                        This method returns the convex relaxation of the hybrid zonotope.
                        If the set is sharp, the convex relaxation is the convex hull.
            
        """
    def copy(self) -> HybZono:
        """copy(self: zonoopt._core.HybZono) -> zonoopt._core.HybZono


                    Creates a copy of the hybrid zonotope object.

                    Returns:
                        HybZono: A copy of the hybrid zonotope object.
            
        """
    def get_A(self) -> scipy.sparse.csc_matrix[numpy.float64]:
        """get_A(self: zonoopt._core.HybZono) -> scipy.sparse.csc_matrix[numpy.float64]


                        Returns constraint matrix
                
                        Returns:
                            scipy.sparse.csc_matrix: A
            
        """
    def get_Ab(self) -> scipy.sparse.csc_matrix[numpy.float64]:
        """get_Ab(self: zonoopt._core.HybZono) -> scipy.sparse.csc_matrix[numpy.float64]


                        Returns binary constraint matrix
                
                        Returns:
                            scipy.sparse.csc_matrix: Ab
            
        """
    def get_Ac(self) -> scipy.sparse.csc_matrix[numpy.float64]:
        """get_Ac(self: zonoopt._core.HybZono) -> scipy.sparse.csc_matrix[numpy.float64]


                        Returns continuous constraint matrix
                
                        Returns:
                            scipy.sparse.csc_matrix: Ac
            
        """
    def get_G(self) -> scipy.sparse.csc_matrix[numpy.float64]:
        """get_G(self: zonoopt._core.HybZono) -> scipy.sparse.csc_matrix[numpy.float64]


                        Returns generator matrix
                
                        Returns:
                            scipy.sparse.csc_matrix: G
            
        """
    def get_Gb(self) -> scipy.sparse.csc_matrix[numpy.float64]:
        """get_Gb(self: zonoopt._core.HybZono) -> scipy.sparse.csc_matrix[numpy.float64]


                        Returns binary generator matrix
                
                        Returns:
                            scipy.sparse.csc_matrix: Gb
            
        """
    def get_Gc(self) -> scipy.sparse.csc_matrix[numpy.float64]:
        """get_Gc(self: zonoopt._core.HybZono) -> scipy.sparse.csc_matrix[numpy.float64]


                        Returns continuous generator matrix
                
                        Returns:
                            scipy.sparse.csc_matrix: Gc
            
        """
    def get_b(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], '[m, 1]']:
        '''get_b(self: zonoopt._core.HybZono) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]


                        Returns constraint vector
                
                        Returns:
                            numpy.array: b
            
        '''
    def get_c(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], '[m, 1]']:
        '''get_c(self: zonoopt._core.HybZono) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]


                        Returns center vector
                
                        Returns:
                            numpy.array: c
            
        '''
    def get_leaves(self, *args, **kwargs):
        """get_leaves(self: zonoopt._core.HybZono, remove_redundancy: bool = True, settings: zonoopt._core.OptSettings = OptSettings structure: verbose: false verbosity_interval: 100 t_max: 1.79769e+308 k_max_admm: 5000 rho: 10 eps_dual: 0.01 eps_prim: 0.001 k_inf_check: 10 inf_norm_conv: true use_interval_contractor: true contractor_iter: 1 search_mode: 0 polish: 1 eps_dual_search: 0.1 eps_prim_search: 0.01 eps_r: 0.01 eps_a: 0.1 k_max_bnb: 100000 n_threads_bnb: 4 n_threads_admm_fp: 3 single_threaded_admm_fp: false max_nodes: 100000 contractor_tree_search_depth: 10 enable_perturb_admm_fp: true k_max_admm_fp_ph1: 10000 k_max_admm_fp_ph2: 90000 cycle_detection_buffer_size: 20 eps_perturb: 0.001 k_restart: 5000 enable_rng_seed: false rng_seed: 0 enable_restart_admm_fp: true, solution: zonoopt._core.OptSolution = None, n_leaves: typing.SupportsInt | typing.SupportsIndex = 2147483647, contractor_iter: typing.SupportsInt | typing.SupportsIndex = 100) -> list[ZonoOpt::ConZono]


                        Computes individual constrained zonotopes whose union is the hybrid zonotope object.
                
                        Args:
                            remove_redundancy (bool, optional): flag to make call to remove_redundancy for each identified leaf
                            settings (OptSettings, optional): optimization settings structure
                            solution (OptSolution, optional): optimization solution structure pointer, populated with result
                            n_leaves (int, optional): max number of leaves to find
                            contractor_iter (int, optional): number of interval contractor iterations to run if using remove_redundancy

                        Returns:
                            list[ConZono]: vector of constrained zonotopes [Z0, Z1, ...] such that Zi is a subset of the current set for all i

                        Searches for constrained zonotopes that correspond to feasible combinations of the hybrid zonotope binary variables.
                        If the branch and bound converges (i.e., did not hit max time, max number of branch and bound iterations, or max nodes in queue)
                        and the n_leaves argument does not stop the optimization before exhausting all possibilities, then the resulting vector of constrained zonotopes
                        can be unioned to recover the original set. It is possible for a leaf to be the empty set if the optimization converges before detecting an infeasibility certificate.
                        Branch and bound search is used to find all leaves of the hybrid zonotope tree. If any threads are allocated
                        for ADMM-FP, these will instead be used for branch and bound search.
            
        """
    def get_n(self) -> int:
        """get_n(self: zonoopt._core.HybZono) -> int


                        Returns dimension of set
                
                        Returns:
                            int: n)
            
        """
    def get_nC(self) -> int:
        """get_nC(self: zonoopt._core.HybZono) -> int


                        Returns number of constraints in set definition
                
                        Returns:
                            int: nC
            
        """
    def get_nG(self) -> int:
        """get_nG(self: zonoopt._core.HybZono) -> int


                        Returns number of generators in set definition
                
                        Returns:
                            int: nG
            
        """
    def get_nGb(self) -> int:
        """get_nGb(self: zonoopt._core.HybZono) -> int


                        Returns number of binary generators in set definition
                
                        Returns:
                            int: nGb
            
        """
    def get_nGc(self) -> int:
        """get_nGc(self: zonoopt._core.HybZono) -> int


                        Returns number of continuous generators in set definition
                
                        Returns:
                            int: nGc
            
        """
    def is_0_1_form(self) -> bool:
        """is_0_1_form(self: zonoopt._core.HybZono) -> bool


                        Returns true if factors are in range [0,1], false if they are in range [-1,1].
                
                        Returns:
                            bool: zero_one_form flag
            
        """
    def is_conzono(self) -> bool:
        """is_conzono(self: zonoopt._core.HybZono) -> bool


                        Polymorphic type checking
                
                        Returns:
                            bool: true if set is a constrained zonotope
            
        """
    def is_empty(self, settings: OptSettings = ..., solution: OptSolution = ..., warm_start_params: WarmStartParams = ...) -> bool:
        """is_empty(self: zonoopt._core.HybZono, settings: zonoopt._core.OptSettings = OptSettings structure: verbose: false verbosity_interval: 100 t_max: 1.79769e+308 k_max_admm: 5000 rho: 10 eps_dual: 0.01 eps_prim: 0.001 k_inf_check: 10 inf_norm_conv: true use_interval_contractor: true contractor_iter: 1 search_mode: 0 polish: 1 eps_dual_search: 0.1 eps_prim_search: 0.01 eps_r: 0.01 eps_a: 0.1 k_max_bnb: 100000 n_threads_bnb: 4 n_threads_admm_fp: 3 single_threaded_admm_fp: false max_nodes: 100000 contractor_tree_search_depth: 10 enable_perturb_admm_fp: true k_max_admm_fp_ph1: 10000 k_max_admm_fp_ph2: 90000 cycle_detection_buffer_size: 20 eps_perturb: 0.001 k_restart: 5000 enable_rng_seed: false rng_seed: 0 enable_restart_admm_fp: true, solution: zonoopt._core.OptSolution = None, warm_start_params: zonoopt._core.WarmStartParams = <zonoopt._core.WarmStartParams object at 0x73c4f34f4570>) -> bool


                        Returns true if the set is provably empty, false otherwise.
                
                        Args:
                            settings (OptSettings, optional): optimization settings structure
                            solution (OptSolution, optional): optimization solution structure pointer, populated with result
                            warm_start_params (WarmStartParams, optional): warm start parameters structure

                        Returns:
                            bool: flag indicating whether set is provably empty
            
        """
    def is_empty_set(self) -> bool:
        """is_empty_set(self: zonoopt._core.HybZono) -> bool


                        Polymorphic type checking

                        Returns:
                            bool: true if set is a empty set object
            
        """
    def is_hybzono(self) -> bool:
        """is_hybzono(self: zonoopt._core.HybZono) -> bool


                        Polymorphic type checking
                
                        Returns:
                            bool: true if set is a hybrid zonotope
            
        """
    def is_point(self) -> bool:
        """is_point(self: zonoopt._core.HybZono) -> bool


                        Polymorphic type checking
                
                        Returns:
                            bool: true if set is a point
            
        """
    def is_sharp(self) -> bool:
        """is_sharp(self: zonoopt._core.HybZono) -> bool


                        Returns true if set is known to be sharp
                
                        Returns:
                            bool: sharp flag
            
        """
    def is_zono(self) -> bool:
        """is_zono(self: zonoopt._core.HybZono) -> bool


                        Polymorphic type checking
                
                        Returns:
                            bool: true if set is a zonotope
            
        """
    def optimize_over(self, P: scipy.sparse.csc_matrix[numpy.float64], q: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], c: typing.SupportsFloat | typing.SupportsIndex = ..., settings: OptSettings = ..., solution: OptSolution = ..., warm_start_params: WarmStartParams = ...) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], '[m, 1]']:
        '''optimize_over(self: zonoopt._core.HybZono, P: scipy.sparse.csc_matrix[numpy.float64], q: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], c: typing.SupportsFloat | typing.SupportsIndex = 0, settings: zonoopt._core.OptSettings = OptSettings structure: verbose: false verbosity_interval: 100 t_max: 1.79769e+308 k_max_admm: 5000 rho: 10 eps_dual: 0.01 eps_prim: 0.001 k_inf_check: 10 inf_norm_conv: true use_interval_contractor: true contractor_iter: 1 search_mode: 0 polish: 1 eps_dual_search: 0.1 eps_prim_search: 0.01 eps_r: 0.01 eps_a: 0.1 k_max_bnb: 100000 n_threads_bnb: 4 n_threads_admm_fp: 3 single_threaded_admm_fp: false max_nodes: 100000 contractor_tree_search_depth: 10 enable_perturb_admm_fp: true k_max_admm_fp_ph1: 10000 k_max_admm_fp_ph2: 90000 cycle_detection_buffer_size: 20 eps_perturb: 0.001 k_restart: 5000 enable_rng_seed: false rng_seed: 0 enable_restart_admm_fp: true, solution: zonoopt._core.OptSolution = None, warm_start_params: zonoopt._core.WarmStartParams = <zonoopt._core.WarmStartParams object at 0x73c4f34f3fb0>) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]


                        Solves optimization problem with quadratic objective over the current set
                
                        Args:
                            P (scipy.sparse.csc_matrix): quadratic objective matrix
                            q (numpy.array): linear objective vector
                            c (float, optional): constant term in objective function
                            settings (OptSettings, optional): optimization settings structure
                            solution (OptSolution, optional): optimization solution structure pointer, populated with result
                            warm_start_params (WarmStartParams, optional): warm start parameters structure

                        Returns:
                            numpy.array: point z in the current set

                        Solves optimization problem of the form min 0.5*z^T*P*z + q^T*z + c where z is a vector in the current set
            
        '''
    def project_point(self, x: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], settings: OptSettings = ..., solution: OptSolution = ..., warm_start_params: WarmStartParams = ...) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], '[m, 1]']:
        '''project_point(self: zonoopt._core.HybZono, x: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], settings: zonoopt._core.OptSettings = OptSettings structure: verbose: false verbosity_interval: 100 t_max: 1.79769e+308 k_max_admm: 5000 rho: 10 eps_dual: 0.01 eps_prim: 0.001 k_inf_check: 10 inf_norm_conv: true use_interval_contractor: true contractor_iter: 1 search_mode: 0 polish: 1 eps_dual_search: 0.1 eps_prim_search: 0.01 eps_r: 0.01 eps_a: 0.1 k_max_bnb: 100000 n_threads_bnb: 4 n_threads_admm_fp: 3 single_threaded_admm_fp: false max_nodes: 100000 contractor_tree_search_depth: 10 enable_perturb_admm_fp: true k_max_admm_fp_ph1: 10000 k_max_admm_fp_ph2: 90000 cycle_detection_buffer_size: 20 eps_perturb: 0.001 k_restart: 5000 enable_rng_seed: false rng_seed: 0 enable_restart_admm_fp: true, solution: zonoopt._core.OptSolution = None, warm_start_params: zonoopt._core.WarmStartParams = <zonoopt._core.WarmStartParams object at 0x73c4f34f4270>) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]


                        Returns the projection of the point x onto the set object.
                
                        Args:
                            x (numpy.array): point to be projected
                            settings (OptSettings, optional): optimization settings structure
                            solution (OptSolution, optional): optimization solution structure pointer, populated with result
                            warm_start_params (WarmStartParams, optional): warm start parameters structure

                        Returns:
                            numpy.array: point z in the current set
            
        '''
    def remove_redundancy(self, contractor_iter: typing.SupportsInt | typing.SupportsIndex = ...) -> bool:
        """remove_redundancy(self: zonoopt._core.HybZono, contractor_iter: typing.SupportsInt | typing.SupportsIndex = 10) -> bool


                        Removes redundant constraints and any unused generators
                
                        This method uses an interval contractor to detect generators that can be removed. 
                        Additionally, any linearly dependent rows of the constraint matrix A are removed.
                        If the linearly dependent constraints are not consistent (e.g., if A = [1, 0.1; 1, 0.1] and b = [1; 0.8]), 
                        the returned set is not equivalent to the original set.
                        Unused factors are also removed.
                
                        Args:
                            contractor_iter (int): number of interval contractor iterations to run

                        Returns:
                            bool: true if successful, false if unable to reduce the complexity of the set representation
            
        """
    def set(self, Gc: scipy.sparse.csc_matrix[numpy.float64], Gb: scipy.sparse.csc_matrix[numpy.float64], c: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], Ac: scipy.sparse.csc_matrix[numpy.float64], Ab: scipy.sparse.csc_matrix[numpy.float64], b: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], zero_one_form: bool = ..., sharp: bool = ...) -> None:
        '''set(self: zonoopt._core.HybZono, Gc: scipy.sparse.csc_matrix[numpy.float64], Gb: scipy.sparse.csc_matrix[numpy.float64], c: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], Ac: scipy.sparse.csc_matrix[numpy.float64], Ab: scipy.sparse.csc_matrix[numpy.float64], b: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], zero_one_form: bool = False, sharp: bool = False) -> None


                        Reset hybrid zonotope object with the given parameters.
                
                        Args:
                            Gc (scipy.sparse.csc_matrix): continuous generator matrix
                            Gb (scipy.sparse.csc_matrix): binary generator matrix
                            c (numpy.array): center
                            Ac (scipy.sparse.csc_matrix): continuous constraint matrix
                            Ab (scipy.sparse.csc_matrix): binary constraint matrix
                            b (numpy.array): constraint vector
                            zero_one_form (bool): true if set is in 0-1 form
                            sharp (bool): true if set is known to be sharp, i.e., convex relaxation = convex hull
            
        '''
    def support(self, d: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], settings: OptSettings = ..., solution: OptSolution = ..., warm_start_params: WarmStartParams = ...) -> float:
        '''support(self: zonoopt._core.HybZono, d: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], settings: zonoopt._core.OptSettings = OptSettings structure: verbose: false verbosity_interval: 100 t_max: 1.79769e+308 k_max_admm: 5000 rho: 10 eps_dual: 0.01 eps_prim: 0.001 k_inf_check: 10 inf_norm_conv: true use_interval_contractor: true contractor_iter: 1 search_mode: 0 polish: 1 eps_dual_search: 0.1 eps_prim_search: 0.01 eps_r: 0.01 eps_a: 0.1 k_max_bnb: 100000 n_threads_bnb: 4 n_threads_admm_fp: 3 single_threaded_admm_fp: false max_nodes: 100000 contractor_tree_search_depth: 10 enable_perturb_admm_fp: true k_max_admm_fp_ph1: 10000 k_max_admm_fp_ph2: 90000 cycle_detection_buffer_size: 20 eps_perturb: 0.001 k_restart: 5000 enable_rng_seed: false rng_seed: 0 enable_restart_admm_fp: true, solution: zonoopt._core.OptSolution = None, warm_start_params: zonoopt._core.WarmStartParams = <zonoopt._core.WarmStartParams object at 0x73c4f34f4770>) -> float


                        Computes support function of the set in the direction d.
                
                        Args:
                            d (numpy.array): vector defining direction for support function
                            settings (OptSettings, optional): optimization settings structure
                            solution (OptSolution, optional): optimization solution structure pointer, populated with result
                            warm_start_params (WarmStartParams, optional): warm start parameters structure

                        Returns:
                            float: support value

                        Solves max_{z in Z} <z, d> where <., .> is the inner product
            
        '''

class IneqTerm:
    """Structure containing term in 0-1 inequality."""
    coeff: float
    idx: int
    def __init__(self, idx: typing.SupportsInt | typing.SupportsIndex, coeff: typing.SupportsFloat | typing.SupportsIndex) -> None:
        """__init__(self: zonoopt._core.IneqTerm, idx: typing.SupportsInt | typing.SupportsIndex, coeff: typing.SupportsFloat | typing.SupportsIndex) -> None

        IneqTerm constructor
        """
    def copy(self) -> IneqTerm:
        """copy(self: zonoopt._core.IneqTerm) -> zonoopt._core.IneqTerm

        Creates a copy of the IneqTerm object.
        """

class IneqType:
    """Enumeration to select inequality direction / use equality.

    Members:

      LESS : Strictly less than

      LESS_OR_EQUAL : Less than or equal to

      EQUAL : Equal to

      GREATER_OR_EQUAL : Greater than or equal to

      GREATER : Strictly greater than"""
    __members__: ClassVar[dict] = ...  # read-only
    EQUAL: ClassVar[IneqType] = ...
    GREATER: ClassVar[IneqType] = ...
    GREATER_OR_EQUAL: ClassVar[IneqType] = ...
    LESS: ClassVar[IneqType] = ...
    LESS_OR_EQUAL: ClassVar[IneqType] = ...
    __entries: ClassVar[dict] = ...
    def __init__(self, value: typing.SupportsInt | typing.SupportsIndex) -> None:
        """__init__(self: zonoopt._core.IneqType, value: typing.SupportsInt | typing.SupportsIndex) -> None"""
    def __eq__(self, other: object) -> bool:
        """__eq__(self: object, other: object, /) -> bool"""
    def __hash__(self) -> int:
        """__hash__(self: object, /) -> int"""
    def __index__(self) -> int:
        """__index__(self: zonoopt._core.IneqType, /) -> int"""
    def __int__(self) -> int:
        """__int__(self: zonoopt._core.IneqType, /) -> int"""
    def __ne__(self, other: object) -> bool:
        """__ne__(self: object, other: object, /) -> bool"""
    @property
    def name(self) -> str:
        """name(self: object, /) -> str

        name(self: object, /) -> str
        """
    @property
    def value(self) -> int:
        """(arg0: zonoopt._core.IneqType) -> int"""

class Inequality:
    """Inequality class"""
    def __init__(self, n_dims: typing.SupportsInt | typing.SupportsIndex) -> None:
        """__init__(self: zonoopt._core.Inequality, n_dims: typing.SupportsInt | typing.SupportsIndex) -> None


                        Constructs inequality, must specify number of dimensions.
                
                        Args:
                            n_dims (int): number of dimensions in the inequality
            
        """
    def add_term(self, idx: typing.SupportsInt | typing.SupportsIndex, coeff: typing.SupportsFloat | typing.SupportsIndex) -> None:
        """add_term(self: zonoopt._core.Inequality, idx: typing.SupportsInt | typing.SupportsIndex, coeff: typing.SupportsFloat | typing.SupportsIndex) -> None


                        Adds a term to the inequality.
                
                        Args:
                            idx (int): index of variable in the inequality
                            coeff (float): coefficient of the variable
            
        """
    def copy(self) -> Inequality:
        """copy(self: zonoopt._core.Inequality) -> zonoopt._core.Inequality


                        Creates a copy of the Inequality object.

                        Returns:
                            Inequality: A copy of the Inequality object.
            
        """
    def get_ineq_type(self) -> IneqType:
        """get_ineq_type(self: zonoopt._core.Inequality) -> zonoopt._core.IneqType


                        Get inequality type / direction.
                
                        Returns:
                            IneqType: inequality type (type member)
            
        """
    def get_n_dims(self) -> int:
        """get_n_dims(self: zonoopt._core.Inequality) -> int


                        Get number of dimensions for inequality.
                
                        Returns:
                            int: number of dimensions (n_dims member)
            
        """
    def get_rhs(self) -> float:
        """get_rhs(self: zonoopt._core.Inequality) -> float


                        Get right-hand side of the inequality.
                
                        Returns:
                            float: right-hand side value (rhs member)
            
        """
    def get_terms(self) -> list[IneqTerm]:
        """get_terms(self: zonoopt._core.Inequality) -> list[zonoopt._core.IneqTerm]


                        Get reference to terms in the inequality.
                
                        Returns:
                            list[IneqTerm]: reference to terms member
            
        """
    def set_ineq_type(self, type: IneqType) -> None:
        """set_ineq_type(self: zonoopt._core.Inequality, type: zonoopt._core.IneqType) -> None


                        Sets the direction of the inequality or sets it to be an equality
                
                        Args:
                            type (IneqType): inequality type (e.g., less than or equal, greater than or equal, or equal)
            
        """
    def set_rhs(self, rhs: typing.SupportsFloat | typing.SupportsIndex) -> None:
        """set_rhs(self: zonoopt._core.Inequality, rhs: typing.SupportsFloat | typing.SupportsIndex) -> None


                        Sets right-hand side of the inequality.
                
                        Args:
                            rhs (float): right-hand side value
            
        """

class Interval:
    """
                Interval class

                Implements interface from IntervalBase. This class owns its lower and upper bounds.
            """
    lb: float
    ub: float
    def __init__(self, y_min: typing.SupportsFloat | typing.SupportsIndex, y_max: typing.SupportsFloat | typing.SupportsIndex) -> None:
        """__init__(self: zonoopt._core.Interval, y_min: typing.SupportsFloat | typing.SupportsIndex, y_max: typing.SupportsFloat | typing.SupportsIndex) -> None


                        Interval constructor.

                        Args:
                            y_min (float): lower bound
                            y_max (float): upper bound
            
        """
    @overload
    def arccos(self) -> Interval:
        """arccos(self: zonoopt._core.Interval) -> zonoopt._core.Interval


                        Compute interval containing arccos(x) for all x in interval

                        Returns:
                            Interval: interval containing arccos(x)
            
        """
    @overload
    def arccos(self, x) -> Any:
        """arccos(self: zonoopt._core.Interval) -> zonoopt._core.Interval


                        Compute interval containing arccos(x) for all x in interval

                        Returns:
                            Interval: interval containing arccos(x)
            
        """
    @overload
    def arccos(self, x) -> Any:
        """arccos(self: zonoopt._core.Interval) -> zonoopt._core.Interval


                        Compute interval containing arccos(x) for all x in interval

                        Returns:
                            Interval: interval containing arccos(x)
            
        """
    @overload
    def arcsin(self) -> Interval:
        """arcsin(self: zonoopt._core.Interval) -> zonoopt._core.Interval


                        Compute interval containing arcsin(x) for all x in interval

                        Returns:
                            Interval: interval containing arcsin(x)
            
        """
    @overload
    def arcsin(self, x) -> Any:
        """arcsin(self: zonoopt._core.Interval) -> zonoopt._core.Interval


                        Compute interval containing arcsin(x) for all x in interval

                        Returns:
                            Interval: interval containing arcsin(x)
            
        """
    @overload
    def arcsin(self, x) -> Any:
        """arcsin(self: zonoopt._core.Interval) -> zonoopt._core.Interval


                        Compute interval containing arcsin(x) for all x in interval

                        Returns:
                            Interval: interval containing arcsin(x)
            
        """
    @overload
    def arctan(self) -> Interval:
        """arctan(self: zonoopt._core.Interval) -> zonoopt._core.Interval


                        Compute interval containing arctan(x) for all x in interval

                        Returns:
                            Interval: interval containing arctan(x)
            
        """
    @overload
    def arctan(self, x) -> Any:
        """arctan(self: zonoopt._core.Interval) -> zonoopt._core.Interval


                        Compute interval containing arctan(x) for all x in interval

                        Returns:
                            Interval: interval containing arctan(x)
            
        """
    @overload
    def arctan(self, x) -> Any:
        """arctan(self: zonoopt._core.Interval) -> zonoopt._core.Interval


                        Compute interval containing arctan(x) for all x in interval

                        Returns:
                            Interval: interval containing arctan(x)
            
        """
    def center(self) -> float:
        """center(self: zonoopt._core.Interval) -> float


                        Gets center of interval (ub + lb) / 2

                        Returns:
                            float: center of interval
            
        """
    def contains(self, y: typing.SupportsFloat | typing.SupportsIndex) -> bool:
        """contains(self: zonoopt._core.Interval, y: typing.SupportsFloat | typing.SupportsIndex) -> bool


                        Checks whether interval contains a value

                        Args:
                            y (float): scalar value

                        Returns:
                            bool: flag indicating if interval contains y
            
        """
    def copy(self) -> Interval:
        """copy(self: zonoopt._core.Interval) -> zonoopt._core.Interval


                        Copy interval object

                        Returns:
                            Interval: copy of interval
            
        """
    @overload
    def cos(self) -> Interval:
        """cos(self: zonoopt._core.Interval) -> zonoopt._core.Interval


                        Compute interval containing cos(x) for all x in interval

                        Returns:
                            Interval: interval containing cos(x)
            
        """
    @overload
    def cos(self, x) -> Any:
        """cos(self: zonoopt._core.Interval) -> zonoopt._core.Interval


                        Compute interval containing cos(x) for all x in interval

                        Returns:
                            Interval: interval containing cos(x)
            
        """
    @overload
    def cos(self, x) -> Any:
        """cos(self: zonoopt._core.Interval) -> zonoopt._core.Interval


                        Compute interval containing cos(x) for all x in interval

                        Returns:
                            Interval: interval containing cos(x)
            
        """
    @overload
    def exp(self) -> Interval:
        """exp(self: zonoopt._core.Interval) -> zonoopt._core.Interval


                        Compute interval containing exp(x) for all x in interval

                        Returns:
                            Interval: interval containing exp(x)
            
        """
    @overload
    def exp(self, x) -> Any:
        """exp(self: zonoopt._core.Interval) -> zonoopt._core.Interval


                        Compute interval containing exp(x) for all x in interval

                        Returns:
                            Interval: interval containing exp(x)
            
        """
    @overload
    def exp(self, x) -> Any:
        """exp(self: zonoopt._core.Interval) -> zonoopt._core.Interval


                        Compute interval containing exp(x) for all x in interval

                        Returns:
                            Interval: interval containing exp(x)
            
        """
    def intersect(self, other: Interval) -> Interval:
        """intersect(self: zonoopt._core.Interval, other: zonoopt._core.Interval) -> zonoopt._core.Interval


                        Interval intersection

                        Args:
                            other (Interval): rhs interval

                        Returns:
                            Interval: intersection of self and other
            
        """
    def inv(self) -> Interval:
        """inv(self: zonoopt._core.Interval) -> zonoopt._core.Interval


                        Interval inverse

                        Returns:
                            Interval: inverse of self
            
        """
    def is_empty(self) -> bool:
        """is_empty(self: zonoopt._core.Interval) -> bool


                        Checks whether interval is empty

                        Returns:
                            bool: flag indicating whether interval is empty
            
        """
    def is_single_valued(self) -> bool:
        """is_single_valued(self: zonoopt._core.Interval) -> bool


                        Checks whether interval is single-valued (i.e., width is 0 within numerical tolerance)

                        Returns:
                            bool: flag indicating if interval is single-value
            
        """
    def radius(self) -> Interval:
        """radius(self: zonoopt._core.Interval) -> zonoopt._core.Interval


                        Gets radius of interval

                        Returns interval centered at zero with width equal to the width of the original interval

                        Returns:
                            Interval: radius of interval
            
        """
    @overload
    def sin(self) -> Interval:
        """sin(self: zonoopt._core.Interval) -> zonoopt._core.Interval


                        Compute interval containing sin(x) for all x in interval

                        Returns:
                            Interval: interval containing sin(x)
            
        """
    @overload
    def sin(self, x) -> Any:
        """sin(self: zonoopt._core.Interval) -> zonoopt._core.Interval


                        Compute interval containing sin(x) for all x in interval

                        Returns:
                            Interval: interval containing sin(x)
            
        """
    @overload
    def sin(self, x) -> Any:
        """sin(self: zonoopt._core.Interval) -> zonoopt._core.Interval


                        Compute interval containing sin(x) for all x in interval

                        Returns:
                            Interval: interval containing sin(x)
            
        """
    @overload
    def tan(self) -> Interval:
        """tan(self: zonoopt._core.Interval) -> zonoopt._core.Interval


                        Compute interval containing tan(x) for all x in interval

                        Returns:
                            Interval: interval containing tan(x)
            
        """
    @overload
    def tan(self, x) -> Any:
        """tan(self: zonoopt._core.Interval) -> zonoopt._core.Interval


                        Compute interval containing tan(x) for all x in interval

                        Returns:
                            Interval: interval containing tan(x)
            
        """
    @overload
    def tan(self, x) -> Any:
        """tan(self: zonoopt._core.Interval) -> zonoopt._core.Interval


                        Compute interval containing tan(x) for all x in interval

                        Returns:
                            Interval: interval containing tan(x)
            
        """
    def width(self) -> float:
        """width(self: zonoopt._core.Interval) -> float


                        Gets width of interval (ub - lb)

                        Returns:
                            float: width of interval
            
        """
    def __add__(self, other: Interval) -> Interval:
        """__add__(self: zonoopt._core.Interval, other: zonoopt._core.Interval) -> zonoopt._core.Interval


                        Interval addition

                        Args:
                            other (Interval): rhs interval

                        Returns:
                            Interval: self + other
            
        """
    @overload
    def __mul__(self, other: Interval) -> Interval:
        """__mul__(*args, **kwargs)
        Overloaded function.

        1. __mul__(self: zonoopt._core.Interval, other: zonoopt._core.Interval) -> zonoopt._core.Interval


                        Interval multiplication

                        Args:
                            other (Interval): rhs interval

                        Returns:
                            Interval: self * other
            

        2. __mul__(self: zonoopt._core.Interval, alpha: typing.SupportsFloat | typing.SupportsIndex) -> zonoopt._core.Interval


                        Interval multiplication with scalar

                        Args:
                            alpha (float): scalar multiplier

                        Returns:
                            Interval: alpha * self
            
        """
    @overload
    def __mul__(self, alpha: typing.SupportsFloat | typing.SupportsIndex) -> Interval:
        """__mul__(*args, **kwargs)
        Overloaded function.

        1. __mul__(self: zonoopt._core.Interval, other: zonoopt._core.Interval) -> zonoopt._core.Interval


                        Interval multiplication

                        Args:
                            other (Interval): rhs interval

                        Returns:
                            Interval: self * other
            

        2. __mul__(self: zonoopt._core.Interval, alpha: typing.SupportsFloat | typing.SupportsIndex) -> zonoopt._core.Interval


                        Interval multiplication with scalar

                        Args:
                            alpha (float): scalar multiplier

                        Returns:
                            Interval: alpha * self
            
        """
    def __sub__(self, other: Interval) -> Interval:
        """__sub__(self: zonoopt._core.Interval, other: zonoopt._core.Interval) -> zonoopt._core.Interval


                        Interval subtraction

                        Args:
                            other (Interval): rhs interval

                        Returns:
                            Interval: self - other
            
        """
    def __truediv__(self, other: Interval) -> Interval:
        """__truediv__(self: zonoopt._core.Interval, other: zonoopt._core.Interval) -> zonoopt._core.Interval


                        Interval division

                        Args:
                            other (Interval): rhs interval

                        Returns:
                            Interval: self / other
            
        """

class IntervalMatrix:
    """Interval matrix class"""
    def __init__(self, mat_lb: scipy.sparse.csc_matrix[numpy.float64], mat_ub: scipy.sparse.csc_matrix[numpy.float64]) -> None:
        """__init__(self: zonoopt._core.IntervalMatrix, mat_lb: scipy.sparse.csc_matrix[numpy.float64], mat_ub: scipy.sparse.csc_matrix[numpy.float64]) -> None


                        IntervalMatrix constructor

                        Args:
                            mat_lb (scipy.sparse.csc_matrix): matrix of lower bounds
                            mat_ub (scipy.sparse.csc_matrix): matrix of upper bounds
            
        """
    def center(self) -> scipy.sparse.csc_matrix[numpy.float64]:
        """center(self: zonoopt._core.IntervalMatrix) -> scipy.sparse.csc_matrix[numpy.float64]


                        Get center matrix

                        Each element of center matrix is the center of the corresponding interval in the interval matrix

                        Returns:
                            scipy.sparse.csc_matrix: center matrix
            
        """
    def cols(self) -> int:
        """cols(self: zonoopt._core.IntervalMatrix) -> int


                        Get number of columns

                        Returns:
                            int: number of cols
            
        """
    def diam(self) -> scipy.sparse.csc_matrix[numpy.float64]:
        """diam(self: zonoopt._core.IntervalMatrix) -> scipy.sparse.csc_matrix[numpy.float64]


                        Get diameter matrix

                        Each element of the diameter matrix is the width of the corresponding interval in the interval matrix

                        Returns:
                            scipy.sparse.csc_matrix: diameter matrix
            
        """
    def radius(self) -> IntervalMatrix:
        """radius(self: zonoopt._core.IntervalMatrix) -> zonoopt._core.IntervalMatrix


                        Get radius of interval matrix

                        Returns interval matrix with intervals centered at zero with width equal to the width of the original interval matrix

                        Returns:
                            IntervalMatrix: radius of interval matrix
            
        """
    def rows(self) -> int:
        """rows(self: zonoopt._core.IntervalMatrix) -> int


                        Get number of rows

                        Returns:
                            int: number of rows
            
        """
    def width(self) -> float:
        """width(self: zonoopt._core.IntervalMatrix) -> float


                        Get width of interval matrix

                        Specifically, this returns the max width for any interval in the interval matrix

                        Returns:
                            float: width of interval matrix
            
        """
    def __add__(self, other: IntervalMatrix) -> IntervalMatrix:
        """__add__(self: zonoopt._core.IntervalMatrix, other: zonoopt._core.IntervalMatrix) -> zonoopt._core.IntervalMatrix


                        IntervalMatrix addition

                        Args:
                            other (IntervalMatrix): rhs interval matrix

                        Returns:
                            IntervalMatrix: resulting interval matrix
            
        """
    @overload
    def __mul__(self, v: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, '[m, 1]']) -> Box:
        '''__mul__(*args, **kwargs)
        Overloaded function.

        1. __mul__(self: zonoopt._core.IntervalMatrix, v: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"]) -> zonoopt._core.Box


                        IntervalMatrix multiplication with vector

                        Args:
                            v (numpy.array): rhs vector

                        Returns:
                            Box: resulting box
            

        2. __mul__(self: zonoopt._core.IntervalMatrix, box: zonoopt._core.Box) -> zonoopt._core.Box


                        IntervalMatrix multiplication with Box

                        Args:
                            box (Box): rhs box

                        Returns:
                            Box: resulting box
            

        3. __mul__(self: zonoopt._core.IntervalMatrix, A: scipy.sparse.csr_matrix[numpy.float64]) -> zonoopt._core.IntervalMatrix


                        IntervalMatrix multiplication with matrix

                        Args:
                            A (scipy.sparse.csr_matrix): rhs matrix

                        Returns:
                            IntervalMatrix: resulting interval matrix
            

        4. __mul__(self: zonoopt._core.IntervalMatrix, other: zonoopt._core.IntervalMatrix) -> zonoopt._core.IntervalMatrix


                        IntervalMatrix multiplication with another IntervalMatrix

                        Args:
                            other (IntervalMatrix): rhs interval matrix

                        Returns:
                            IntervalMatrix: resulting interval matrix
            
        '''
    @overload
    def __mul__(self, box: Box) -> Box:
        '''__mul__(*args, **kwargs)
        Overloaded function.

        1. __mul__(self: zonoopt._core.IntervalMatrix, v: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"]) -> zonoopt._core.Box


                        IntervalMatrix multiplication with vector

                        Args:
                            v (numpy.array): rhs vector

                        Returns:
                            Box: resulting box
            

        2. __mul__(self: zonoopt._core.IntervalMatrix, box: zonoopt._core.Box) -> zonoopt._core.Box


                        IntervalMatrix multiplication with Box

                        Args:
                            box (Box): rhs box

                        Returns:
                            Box: resulting box
            

        3. __mul__(self: zonoopt._core.IntervalMatrix, A: scipy.sparse.csr_matrix[numpy.float64]) -> zonoopt._core.IntervalMatrix


                        IntervalMatrix multiplication with matrix

                        Args:
                            A (scipy.sparse.csr_matrix): rhs matrix

                        Returns:
                            IntervalMatrix: resulting interval matrix
            

        4. __mul__(self: zonoopt._core.IntervalMatrix, other: zonoopt._core.IntervalMatrix) -> zonoopt._core.IntervalMatrix


                        IntervalMatrix multiplication with another IntervalMatrix

                        Args:
                            other (IntervalMatrix): rhs interval matrix

                        Returns:
                            IntervalMatrix: resulting interval matrix
            
        '''
    @overload
    def __mul__(self, A: scipy.sparse.csr_matrix[numpy.float64]) -> IntervalMatrix:
        '''__mul__(*args, **kwargs)
        Overloaded function.

        1. __mul__(self: zonoopt._core.IntervalMatrix, v: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"]) -> zonoopt._core.Box


                        IntervalMatrix multiplication with vector

                        Args:
                            v (numpy.array): rhs vector

                        Returns:
                            Box: resulting box
            

        2. __mul__(self: zonoopt._core.IntervalMatrix, box: zonoopt._core.Box) -> zonoopt._core.Box


                        IntervalMatrix multiplication with Box

                        Args:
                            box (Box): rhs box

                        Returns:
                            Box: resulting box
            

        3. __mul__(self: zonoopt._core.IntervalMatrix, A: scipy.sparse.csr_matrix[numpy.float64]) -> zonoopt._core.IntervalMatrix


                        IntervalMatrix multiplication with matrix

                        Args:
                            A (scipy.sparse.csr_matrix): rhs matrix

                        Returns:
                            IntervalMatrix: resulting interval matrix
            

        4. __mul__(self: zonoopt._core.IntervalMatrix, other: zonoopt._core.IntervalMatrix) -> zonoopt._core.IntervalMatrix


                        IntervalMatrix multiplication with another IntervalMatrix

                        Args:
                            other (IntervalMatrix): rhs interval matrix

                        Returns:
                            IntervalMatrix: resulting interval matrix
            
        '''
    @overload
    def __mul__(self, other: IntervalMatrix) -> IntervalMatrix:
        '''__mul__(*args, **kwargs)
        Overloaded function.

        1. __mul__(self: zonoopt._core.IntervalMatrix, v: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"]) -> zonoopt._core.Box


                        IntervalMatrix multiplication with vector

                        Args:
                            v (numpy.array): rhs vector

                        Returns:
                            Box: resulting box
            

        2. __mul__(self: zonoopt._core.IntervalMatrix, box: zonoopt._core.Box) -> zonoopt._core.Box


                        IntervalMatrix multiplication with Box

                        Args:
                            box (Box): rhs box

                        Returns:
                            Box: resulting box
            

        3. __mul__(self: zonoopt._core.IntervalMatrix, A: scipy.sparse.csr_matrix[numpy.float64]) -> zonoopt._core.IntervalMatrix


                        IntervalMatrix multiplication with matrix

                        Args:
                            A (scipy.sparse.csr_matrix): rhs matrix

                        Returns:
                            IntervalMatrix: resulting interval matrix
            

        4. __mul__(self: zonoopt._core.IntervalMatrix, other: zonoopt._core.IntervalMatrix) -> zonoopt._core.IntervalMatrix


                        IntervalMatrix multiplication with another IntervalMatrix

                        Args:
                            other (IntervalMatrix): rhs interval matrix

                        Returns:
                            IntervalMatrix: resulting interval matrix
            
        '''
    def __sub__(self, other: IntervalMatrix) -> IntervalMatrix:
        """__sub__(self: zonoopt._core.IntervalMatrix, other: zonoopt._core.IntervalMatrix) -> zonoopt._core.IntervalMatrix


                        IntervalMatrix subtraction

                        Args:
                            other (IntervalMatrix): rhs interval matrix

                        Returns:
                            IntervalMatrix: resulting interval matrix
            
        """

class OptSettings:
    """Settings for optimization routines in ZonoOpt library."""
    contractor_iter: int
    contractor_tree_search_depth: int
    cycle_detection_buffer_size: int
    enable_perturb_admm_fp: bool
    enable_restart_admm_fp: bool
    enable_rng_seed: bool
    eps_a: float
    eps_dual: float
    eps_dual_search: float
    eps_perturb: float
    eps_prim: float
    eps_prim_search: float
    eps_r: float
    inf_norm_conv: bool
    k_inf_check: int
    k_max_admm: int
    k_max_admm_fp_ph1: int
    k_max_admm_fp_ph2: int
    k_max_bnb: int
    k_restart: int
    max_nodes: int
    n_threads_admm_fp: int
    n_threads_bnb: int
    polish: bool
    rho: float
    rng_seed: int
    search_mode: int
    single_threaded_admm_fp: bool
    t_max: float
    use_interval_contractor: bool
    verbose: bool
    verbosity_interval: int
    def __init__(self) -> None:
        """__init__(self: zonoopt._core.OptSettings) -> None"""
    def copy(self) -> OptSettings:
        """copy(self: zonoopt._core.OptSettings) -> zonoopt._core.OptSettings


                        Copy settings object

                        Returns:
                            OptSettings: copy of settings
            
        """
    def settings_valid(self) -> bool:
        """settings_valid(self: zonoopt._core.OptSettings) -> bool

        check whether settings struct is valid
        """

class OptSolution:
    """Solution data structure for optimization routines in ZonoOpt library."""
    J: float
    converged: bool
    dual_residual: float
    infeasible: bool
    iter: int
    primal_residual: float
    run_time: float
    startup_time: float
    u: typing.Annotated[numpy.typing.NDArray[numpy.float64], '[m, 1]']
    x: typing.Annotated[numpy.typing.NDArray[numpy.float64], '[m, 1]']
    z: typing.Annotated[numpy.typing.NDArray[numpy.float64], '[m, 1]']
    def __init__(self) -> None:
        """__init__(self: zonoopt._core.OptSolution) -> None"""
    def copy(self) -> OptSolution:
        """copy(self: zonoopt._core.OptSolution) -> zonoopt._core.OptSolution


                        Copy solution object

                        Returns:
                            OptSolution: copy of solution
            
        """

class Point(Zono):
    """
                Point class
                
                A point is defined entirely by the center vector c.
            """
    def __init__(self, c: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, '[m, 1]']) -> None:
        '''__init__(self: zonoopt._core.Point, c: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"]) -> None


                        Point constructor
                
                        Args:
                            c (numpy.array): center vector
            
        '''
    def set(self, c: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, '[m, 1]']) -> None:
        '''set(self: zonoopt._core.Point, c: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"]) -> None


                        Reset point object with the given parameters.
                
                        Args:
                            c (numpy.array): center vector
            
        '''

class WarmStartParams:
    """
            Warm start parameters for optimization routines in ZonoOpt library.
        """
    u: typing.Annotated[numpy.typing.NDArray[numpy.float64], '[m, 1]']
    z: typing.Annotated[numpy.typing.NDArray[numpy.float64], '[m, 1]']
    def __init__(self) -> None:
        """__init__(self: zonoopt._core.WarmStartParams) -> None"""
    def copy(self) -> WarmStartParams:
        """copy(self: zonoopt._core.WarmStartParams) -> zonoopt._core.WarmStartParams


                        Copy warm start parameters object

                        Returns:
                            WarmStartParams: copy of warm start parameters
            
        """

class Zono(ConZono):
    """
                Zonotope class
                
                A zonotope is defined as:
                Z = {G \\xi + c | \\xi in [-1, 1]^nG}.
                Equivalently, the following shorthand can be used: Z = <G, c>.
                Optionally, in 0-1 form, the factors are xi in [0,1].
                The set dimension is n, and the number of generators is nG.
            """
    def __init__(self, G: scipy.sparse.csc_matrix[numpy.float64], c: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], zero_one_form: bool = ...) -> None:
        '''__init__(self: zonoopt._core.Zono, G: scipy.sparse.csc_matrix[numpy.float64], c: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], zero_one_form: bool = False) -> None


                        Zono constructor
                
                        Args:
                            G (scipy.sparse.csc_matrix): generator matrix
                            c (numpy.array): center
                            zero_one_form (bool, optional): true if set is in 0-1 form
            
        '''
    def get_center(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], '[m, 1]']:
        '''get_center(self: zonoopt._core.Zono) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]


                        Get center of zonotope.

                        Returns:
                            numpy.array: center vector
            
        '''
    def get_volume(self) -> float:
        '''get_volume(self: zonoopt._core.Zono) -> float


                            Get volume of zonotope.

                            Reference: Gover and Krikorian 2010, "Determinants and the volumes of parallelotopes and zonotopes"
                            Requires nG choose n determinant computations.

                            Returns:
                                float: volume of zonotope
            
        '''
    def reduce_order(self, n_o: typing.SupportsInt | typing.SupportsIndex) -> Zono:
        """reduce_order(self: zonoopt._core.Zono, n_o: typing.SupportsInt | typing.SupportsIndex) -> zonoopt._core.Zono


                        Perform zonotope order reduction.

                        Args:
                            n_o (int): desired order, must be greater than or equal to the dimension of the set

                        Returns:
                            zonotope with order n_o
            
        """
    def set(self, G: scipy.sparse.csc_matrix[numpy.float64], c: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], zero_one_form: bool = ...) -> None:
        '''set(self: zonoopt._core.Zono, G: scipy.sparse.csc_matrix[numpy.float64], c: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], zero_one_form: bool = False) -> None


                        Reset zonotope object with the given parameters.
                
                        Args:
                            G (scipy.sparse.csc_matrix): generator matrix
                            c (numpy.array): center
                            zero_one_form (bool, optional): true if set is in 0-1 form
            
        '''

def affine_inclusion(Z: HybZono, R: IntervalMatrix, s: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'] = ...) -> HybZono:
    '''affine_inclusion(Z: zonoopt._core.HybZono, R: zonoopt._core.IntervalMatrix, s: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"] = array([], dtype=float64)) -> zonoopt._core.HybZono


                Returns inclusion of zonotopic set for uncertain affine map R*Z + s

                This computes an over-approximation of the affine map using the method of
                Rego et. al. (2020) "Guaranteed methods based on constrained zonotopes for set-valued state estimation of nonlinear discrete-time systems"
                The SVD-based zonotope over-approximation method is used in this function when Z is a constrained zonotope.
                When Z is a hybrid zonotope, the convex relaxation is used to produce a constrained zonotope, and then the SVD-based method is applied.

                Args:
                    Z (HybZono): zonotopic set
                    R (IntervalMatrix): affine map interval matrix
                    s (numpy.array, optional): vector offset

                Returns:
                    HybZono: zonotopic set
        
    '''
def affine_map(Z: HybZono, R: scipy.sparse.csc_matrix[numpy.float64], s: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'] = ...) -> HybZono:
    '''affine_map(Z: zonoopt._core.HybZono, R: scipy.sparse.csc_matrix[numpy.float64], s: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"] = array([], dtype=float64)) -> zonoopt._core.HybZono


                Returns affine map R*Z + s of set Z
            
                Args:
                    Z (HybZono): zonotopic set
                    R (scipy.sparse.csc_matrix): affine map matrix
                    s (numpy.array, optional): vector offset
            
                Returns:
                    HybZono: zonotopic set
        
    '''
def cartesian_product(Z1: HybZono, Z2: HybZono) -> HybZono:
    """cartesian_product(Z1: zonoopt._core.HybZono, Z2: zonoopt._core.HybZono) -> zonoopt._core.HybZono


                Computes the Cartesian product of two sets Z1 and Z2.
            
                Args:
                    Z1 (HybZono): zonotopic set
                    Z2 (HybZono): zonotopic set
            
                Returns:
                    HybZono: zonotopic set
        
    """
def constrain(Z: HybZono, ineqs: collections.abc.Sequence[Inequality], R: scipy.sparse.csc_matrix[numpy.float64] = ...) -> HybZono:
    """constrain(Z: zonoopt._core.HybZono, ineqs: collections.abc.Sequence[zonoopt._core.Inequality], R: scipy.sparse.csc_matrix[numpy.float64] = <Compressed Sparse Column sparse matrix of dtype 'float64' with 0 stored elements and shape (0, 0)>) -> zonoopt._core.HybZono


                Applies inequalities to set.
            
                Args:
                    Z (HybZono): Set for inequalities to be applied to
                    ineqs (list[Inequality]): list of inequalities
                    R (scipy.sparse.csc_matrix, optional): For generalized intersection with the inequalities. Defaults to identity.
            
                Returns:
                    HybZono: zonotopic set

                Constrains set Z by applying the given inequalities to the set.
                For example, given a set Z with states z0, z1, z2, the constraint z0 + z1 - z2 <= 2 could be added via an inequality object.
                R is used for generalized intersection-like operations. For instance, when all the inequalities are <= inequalities,
                this function returns Z int_R (Hx<=f) where H is the halfspace represented by the inequalities.
        
    """
def convex_hull(Zs: collections.abc.Sequence[HybZono]) -> ConZono:
    """convex_hull(Zs: collections.abc.Sequence[zonoopt._core.HybZono]) -> zonoopt._core.ConZono


                Computes the convex hull of several sets

                Computes convex hull of sets {Z0, Z1, ..., Zn}.
                If Zi is a hybrid zonotope, it must be sharp or this function will throw an error.

                Args:
                    Zs (list[HybZono]): sets for which convex hull is to be computed.

                Returns:
                    ConZono: constrained zonotop convex hull
        
    """
def halfspace_intersection(Z: HybZono, H: scipy.sparse.csc_matrix[numpy.float64], f: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], R: scipy.sparse.csc_matrix[numpy.float64] = ...) -> HybZono:
    '''halfspace_intersection(Z: zonoopt._core.HybZono, H: scipy.sparse.csc_matrix[numpy.float64], f: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], R: scipy.sparse.csc_matrix[numpy.float64] = <Compressed Sparse Column sparse matrix of dtype \'float64\' with 0 stored elements and shape (0, 0)>) -> zonoopt._core.HybZono


                Computes the intersection generalized intersection of set Z with halfspace H*x <= f over matrix R.
            
                Args:
                    Z (HybZono): zonotopic set
                    H (scipy.sparse.csc_matrix): halfspace matrix
                    f (numpy.array): halfspace vector
                    R (scipy.sparse.csc_matrix, optional): affine map matrix
            
                Returns:
                    HybZono: zonotopic set
        
    '''
def intersection(Z1: HybZono, Z2: HybZono, R: scipy.sparse.csc_matrix[numpy.float64] = ...) -> HybZono:
    """intersection(Z1: zonoopt._core.HybZono, Z2: zonoopt._core.HybZono, R: scipy.sparse.csc_matrix[numpy.float64] = <Compressed Sparse Column sparse matrix of dtype 'float64' with 0 stored elements and shape (0, 0)>) -> zonoopt._core.HybZono


                Computes the generalized intersection of sets Z1 and Z2 over the matrix R.
            
                Args:
                    Z1 (HybZono): zonotopic set
                    Z2 (HybZono): zonotopic set
                    R (scipy.sparse.csc_matrix, optional): affine map matrix
            
                Returns:
                    HybZono: zonotopic set
        
    """
def intersection_over_dims(Z1: HybZono, Z2: HybZono, dims: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> HybZono:
    """intersection_over_dims(Z1: zonoopt._core.HybZono, Z2: zonoopt._core.HybZono, dims: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> zonoopt._core.HybZono


                Computes the intersection of sets Z1 and Z2 over the specified dimensions.
            
                Args:
                    Z1 (HybZono): zonotopic set
                    Z2 (HybZono): zonotopic set
                    dims (list[int]): list of dimensions
            
                Returns:
                    HybZono: zonotopic set
        
    """
def interval_2_zono(box: Box) -> Zono:
    """interval_2_zono(box: zonoopt._core.Box) -> zonoopt._core.Zono


                Builds a zonotope from a Box object.
            
                Args:
                    box (Box): Box object (vector of intervals)
            
                Returns:
                    Zono: zonotope
        
    """
def make_regular_zono_2D(radius: typing.SupportsFloat | typing.SupportsIndex, n_sides: typing.SupportsInt | typing.SupportsIndex, outer_approx: bool = ..., c: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, '[2, 1]'] = ...) -> Zono:
    '''make_regular_zono_2D(radius: typing.SupportsFloat | typing.SupportsIndex, n_sides: typing.SupportsInt | typing.SupportsIndex, outer_approx: bool = False, c: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[2, 1]"] = array([0., 0.])) -> zonoopt._core.Zono


                Builds a 2D regular zonotope with a given radius and number of sides.
            
                Args:
                    radius (float): radius of the zonotope
                    n_sides (int): number of sides (must be an even number >= 4)
                    outer_approx (bool, optional): flag to do an outer approximation instead of an inner approximation
                    c (numpy.array, optional): center vector
            
                Returns:
                    Zono: zonotope
        
    '''
def minkowski_sum(Z1: HybZono, Z2: HybZono) -> HybZono:
    """minkowski_sum(Z1: zonoopt._core.HybZono, Z2: zonoopt._core.HybZono) -> zonoopt._core.HybZono


                Computes Minkowski sum of two sets Z1 and Z2.
            
                Args:
                    Z1 (HybZono): zonotopic set
                    Z2 (HybZono): zonotopic set
            
                Returns:
                    HybZono: zonotopic set
        
    """
def pontry_diff(Z1: HybZono, Z2: HybZono, exact: bool = ...) -> HybZono:
    """pontry_diff(Z1: zonoopt._core.HybZono, Z2: zonoopt._core.HybZono, exact: bool = False) -> zonoopt._core.HybZono


                Computes the Pontryagin difference Z1 - Z2.
            
                For inner approximations (exact=false), the algorithm from Vinod et. al. 2025 is used.
                Note that this algorithm is exact when the minuend is a constrained zonotope and the matrix [G;A] is invertible.
                Exact Pontryagin difference can only be computed when the subtrahend is a zonotope.
                If subtrahend is a constrained zonotope, it will first be over-approximated as a zonotope.
                If subtrahend is a hybrid zonotope, a get_leaves operation will first be performed to produce
                a union of constrained zonotopes.
                If the minuend is a hybrid zonotope and exact is false, a get_leaves operation will be performed for 
                Z1 to reduce the number of leaves in the resulting set.
            
                Args:
                    Z1 (HybZono): minuend
                    Z2 (ConZono): subtrahend
                    exact (bool, optional): require output to be exact, otherwise inner approximation will be returned
            
                Returns:
                    HybZono: zonotopic set
        
    """
def project_onto_dims(Z: HybZono, dims: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> HybZono:
    """project_onto_dims(Z: zonoopt._core.HybZono, dims: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> zonoopt._core.HybZono


                Projects set Z onto the dimensions specified in dims.
            
                Args:
                    Z (HybZono): zonotopic set
                    dims (list[int]): list of dimensions to project onto
            
                Returns:
                    HybZono: zonotopic set
        
    """
def set_diff(Z1: HybZono, Z2: HybZono, delta_m: typing.SupportsFloat | typing.SupportsIndex = ..., remove_redundancy: bool = ..., settings: OptSettings = ..., solution: OptSolution = ..., n_leaves: typing.SupportsInt | typing.SupportsIndex = ..., contractor_iter: typing.SupportsInt | typing.SupportsIndex = ...) -> HybZono:
    """set_diff(Z1: zonoopt._core.HybZono, Z2: zonoopt._core.HybZono, delta_m: typing.SupportsFloat | typing.SupportsIndex = 100, remove_redundancy: bool = True, settings: zonoopt._core.OptSettings = OptSettings structure: verbose: false verbosity_interval: 100 t_max: 1.79769e+308 k_max_admm: 5000 rho: 10 eps_dual: 0.01 eps_prim: 0.001 k_inf_check: 10 inf_norm_conv: true use_interval_contractor: true contractor_iter: 1 search_mode: 0 polish: 1 eps_dual_search: 0.1 eps_prim_search: 0.01 eps_r: 0.01 eps_a: 0.1 k_max_bnb: 100000 n_threads_bnb: 4 n_threads_admm_fp: 3 single_threaded_admm_fp: false max_nodes: 100000 contractor_tree_search_depth: 10 enable_perturb_admm_fp: true k_max_admm_fp_ph1: 10000 k_max_admm_fp_ph2: 90000 cycle_detection_buffer_size: 20 eps_perturb: 0.001 k_restart: 5000 enable_rng_seed: false rng_seed: 0 enable_restart_admm_fp: true, solution: zonoopt._core.OptSolution = None, n_leaves: typing.SupportsInt | typing.SupportsIndex = 2147483647, contractor_iter: typing.SupportsInt | typing.SupportsIndex = 10) -> zonoopt._core.HybZono


                Set difference Z1 \\\\ Z2

                Args:
                    Z1 (HybZono): zonotopic set
                    Z2 (HybZono): zonotopic set
                    delta_m (float, optional): parameter defining range of complement
                    remove_redundancy (bool, optional): remove redundant constraints and unused generators in get_leaves function call
                    settings (OptSettings, optional): optimization settings for get_leaves function call
                    solution (OptSolution, optional): optimization solution for get_leaves function call
                    n_leaves (int, optional): maximum number of leaves to return in get_leaves function call
                    contractor_iter (int, optional): number of interval contractor iterations if using remove_redundancy

                Returns:
                    HybZono: zonotopic set
        
    """
def union_of_many(Z_list: collections.abc.Sequence[HybZono], preserve_sharpness: bool = ..., expose_indicators: bool = ...) -> HybZono:
    """union_of_many(Z_list: collections.abc.Sequence[zonoopt._core.HybZono], preserve_sharpness: bool = False, expose_indicators: bool = False) -> zonoopt._core.HybZono


                Computes the union of several sets
            
                Args:
                    Z_list (list[HybZono]): sets to be unioned
                    preserve_sharpness (bool, optional): flag to preserve sharpness of the union at expense of complexity.
                    expose_indicators (bool, optional): flag to append indicator set to the union.
            
                Returns:
                    HybZono: zonotopic set
        
                Computes union of sets {Z0, Z1, ..., Zn}. If expose_indicators is true, returns union({Z0, ..., Zn}) x I where I is the indicator set for the union.
                Specifically, each dimension of I corresponds to one of the Zi in the union. So for union_of_many({Z0, Z1, Z2}, true) with Z0, Z1, Z2 not intersecting,
                if a vector [z, i] is in union({Z0, Z1, Z2}) x I, then i = [1, 0, 0] if z is in Z0, etc.
        
    """
def vrep_2_conzono(Vpoly: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, '[m, n]']) -> ConZono:
    '''vrep_2_conzono(Vpoly: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, n]"]) -> zonoopt._core.ConZono


                Builds a constrained zonotope from a vertex representation polytope.

                Args:
                    Vpoly (numpy.array): vertices of V-rep polytope
            
                Returns:
                    ConZono: constrained zonotope
            
                Vpoly is a matrix where each row is a vertex of the polytope.
        
    '''
def vrep_2_hybzono(Vpolys: collections.abc.Sequence[typing.Annotated[numpy.typing.ArrayLike, numpy.float64, '[m, n]']], expose_indicators: bool = ...) -> HybZono:
    '''vrep_2_hybzono(Vpolys: collections.abc.Sequence[typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, n]"]], expose_indicators: bool = False) -> zonoopt._core.HybZono


                Computes a hybrid zonotope from a union of vertex representation polytopes.

                Args:
                    Vpolys (list[numpy.array]): V-rep polytopes to be unioned
                    expose_indicators (bool, optional): flag to append indicator set to the union.
            
                Returns:
                    HybZono: hybrid zonotope
            
                Vpolys is a vector of matrices, where each matrix represents a polytope in vertex representation.
                Each row in each polytope matrix is a vertex of the polytope, and each column corresponds to a dimension.
                The function constructs a hybrid zonotope in [0,1] form that represents the union of these polytopes.
                This function computes union of sets {V0, V1, ..., Vn}. If expose_indicators is true, returns union({V0, ..., Vn}) x I where I is the indicator set for the union.
                Specifically, each dimension of I corresponds to one of the Vi in the union. So for vrep_2_hybzono({V0, V1, V2}, true) with V0, V1, V2 not intersecting,
                if a vector [z, i] is in union({V0, V1, V2}) x I, then i = [1, 0, 0] if z is in V0, etc.
        
    '''
def zono_union_2_hybzono(Zs: collections.abc.Sequence[Zono], expose_indicators: bool = ...) -> HybZono:
    """zono_union_2_hybzono(Zs: collections.abc.Sequence[zonoopt._core.Zono], expose_indicators: bool = False) -> zonoopt._core.HybZono


                Computes a hybrid zonotope from a union of zonotopes.

                Args:
                    Zs (list[Zono]): zonotopes to be unioned
                    expose_indicators (bool, optional): flag to append indicator set to the union.
            
                Returns:
                    HybZono: hybrid zonotope
            
                This function computes union of sets {Z0, Z1, ..., Zn}. This can be more efficient than union_of_many if all sets are zonotopes because generators can be reused.
                If expose_indicators is true, returns union({Z0, ..., Zn}) x I where I is the indicator set for the union.
                Specifically, each dimension of I corresponds to one of the Zi in the union. So for zono_union_2_hybzono({Z0, Z1, Z2}, true) with Z0, Z1, VZ2 not intersecting,
                if a vector [z, i] is in union({Z0, Z1, Z2}) x I, then i = [1, 0, 0] if z is in Z0, etc.
        
    """
