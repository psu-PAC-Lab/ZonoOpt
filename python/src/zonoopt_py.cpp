#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#define IS_PYTHON_ENV
#define zono_float double
#include "ZonoOpt.hpp"
using namespace ZonoOpt;

#include <sstream>

PYBIND11_MODULE(_core, m)
{
    m.doc() = R"pbdoc(
        Classes and tailored optimization routines for zonotopes, constrained zonotopes, and hybrid zonotopes.

        See the github page for documentation: https://github.com/psu-PAC-Lab/ZonoOpt

        More information about ZonoOpt can be found in the the following publication. Please cite this if you publish work based on ZonoOpt: 
        Robbins, J.A., Siefert, J.A., and Pangborn, H.C., "Sparsity-Promoting Reachability Analysis and Optimization of Constrained Zonotopes," 2025.**
    )pbdoc";

    // solver settings and solution classes
    py::class_<OptSettings>(m, "OptSettings", "Settings for optimization routines in ZonoOpt library.")
        .def(py::init())
        .def_readwrite("verbose", &OptSettings::verbose, "display optimization progress")
        .def_readwrite("verbosity_interval", &OptSettings::verbosity_interval, "print every verbose_interval iterations")
        .def_readwrite("t_max", &OptSettings::t_max, "max time for optimization")
        .def_readwrite("use_interval_contractor", &OptSettings::use_interval_contractor, "flag to use interval contractor for constraint tightening / implication")
        .def_readwrite("contractor_iter", &OptSettings::contractor_iter, "number of interval contractor iterations")
        .def_readwrite("k_max_admm", &OptSettings::k_max_admm, "max convex admm iterations")
        .def_readwrite("rho", &OptSettings::rho, "admm penalty parameter, higher prioritizes feasibility during iterations, lower prioritizes optimality")
        .def_readwrite("eps_dual", &OptSettings::eps_dual, "dual convergence tolerance")
        .def_readwrite("eps_prim", &OptSettings::eps_prim, "primal convergence tolerance")
        .def_readwrite("k_inf_check", &OptSettings::k_inf_check, "check infeasibility every k_inf_check iterations")
        .def_readwrite("inf_norm_conv", &OptSettings::inf_norm_conv, 
            "use infinity norm for convergence check (if false, scaled 2-norm is used)")
        .def_readwrite("search_mode", &OptSettings::search_mode, 
            "0-> best first, 1-> best dive")
        .def_readwrite("polish", &OptSettings::polish, "flag to perform solution polishing")
        .def_readwrite("eps_dual_search", &OptSettings::eps_dual_search, 
            "dual residual convergence tolerance during branch and bound and search")
        .def_readwrite("eps_prim_search", &OptSettings::eps_prim_search,
            "primal residual convergence tolerance during branch and bound and search")
        .def_readwrite("eps_r", &OptSettings::eps_r, "relative convergence tolerance")
        .def_readwrite("eps_a", &OptSettings::eps_a, "absolute convergence tolerance")
        .def_readwrite("k_max_bnb", &OptSettings::k_max_bnb, "max number of branch-and-bound iterations")
        .def_readwrite("n_threads_bnb", &OptSettings::n_threads_bnb, "max threads for branch and bound")
        .def_readwrite("n_threads_admm_fp", &OptSettings::n_threads_admm_fp, "max threads for ADMM-FP")
        .def_readwrite("single_threaded_admm_fp", &OptSettings::single_threaded_admm_fp,
            "enables single-threaded ADMM-FP solution, overrides n_threads_bnb, n_threads_admm_fp")
        .def_readwrite("max_nodes", &OptSettings::max_nodes, "terminate if more than this many nodes are in branch and bound queue")
        .def_readwrite("contractor_tree_search_depth", &OptSettings::contractor_tree_search_depth, 
            "when applying interval contractor in branch and bound, this is how deep to search the constraint tree for affected variables")
        .def_readwrite("k_max_admm_fp_ph1", &OptSettings::k_max_admm_fp_ph1,
            "max ADMM iterations for ADMM-FP phase 1 (objective included)")
        .def_readwrite("k_max_admm_fp_ph2", &OptSettings::k_max_admm_fp_ph2,
            "max ADMM iterations for ADMM-FP phase 2 (no objective)")
        .def_readwrite("cycle_detection_buffer_size", &OptSettings::cycle_detection_buffer_size,
            "in ADMM-FP, this is the max size of the buffer that checks for cycles")
        .def_readwrite("enable_perturb_admm_fp", &OptSettings::enable_perturb_admm_fp,
            "enable perturbations in ADMM-FP")
        .def_readwrite("enable_restart_admm_fp", &OptSettings::enable_restart_admm_fp,
            "enable restarts (significant perturbations) in ADMM-FP")
        .def_readwrite("eps_perturb", &OptSettings::eps_perturb,
            "relative tolerance for cycling detection in ADMM-FP, triggers perturbation")
        .def_readwrite("k_restart", &OptSettings::k_restart,
            "perform restart operation if primal residual does not improve over this many iterations in ADMM-FP")
        .def_readwrite("enable_rng_seed", &OptSettings::enable_rng_seed,
            "enable rng seed for ADMM-FP")
        .def_readwrite("rng_seed", &OptSettings::rng_seed,
            "rng seed for ADMM-FP")
        .def("settings_valid", &OptSettings::settings_valid, "check whether settings struct is valid")
        .def("__repr__", &OptSettings::print,
            R"pbdoc(
                Displays settings as a string

                Returns:
                    str: string
            )pbdoc")
        .def("copy", [](const  OptSettings& self) -> OptSettings
            { return self; },
            R"pbdoc(
                Copy settings object

                Returns:
                    OptSettings: copy of settings
            )pbdoc")
    ;

    py::class_<OptSolution>(m, "OptSolution", "Solution data structure for optimization routines in ZonoOpt library.")
        .def(py::init())
        .def_readwrite("z", &OptSolution::z, "solution vector")
        .def_readwrite("J", &OptSolution::J, "objective")
        .def_readwrite("run_time", &OptSolution::run_time, "time to compute solution")
        .def_readwrite("startup_time", &OptSolution::startup_time, "time to factorize matrices and run interval contractors")
        .def_readwrite("iter", &OptSolution::iter, "number of iterations")
        .def_readwrite("converged", &OptSolution::converged, "true if optimization has converged")
        .def_readwrite("infeasible", &OptSolution::infeasible, "true if optimization problem is provably infeasible")
        .def_readwrite("x", &OptSolution::x, "ADMM primal variable, approximately equal to z when converged")
        .def_readwrite("u", &OptSolution::u, "ADMM dual variable")
        .def_readwrite("primal_residual", &OptSolution::primal_residual, "primal residual, corresponds to feasibility")
        .def_readwrite("dual_residual", &OptSolution::dual_residual, "dual residual, corresponds to optimality")
        .def("__repr__", &OptSolution::print,
            R"pbdoc(
                Displays solution as a string

                Returns:
                    str: string
            )pbdoc")
        .def("copy", [](const OptSolution& self) -> OptSolution
            { return self; },
            R"pbdoc(
                Copy solution object

                Returns:
                    OptSolution: copy of solution
            )pbdoc")
    ;

    py::class_<WarmStartParams>(m, "WarmStartParams",
        R"pbdoc(
            Warm start parameters for optimization routines in ZonoOpt library.
        )pbdoc")
        .def(py::init())
        .def_readwrite("z", &WarmStartParams::z, "warm start primal variable")
        .def_readwrite("u", &WarmStartParams::u, "warm start dual variable")
        .def("copy", [](const WarmStartParams& self) -> WarmStartParams
            { return self; },
            R"pbdoc(
                Copy warm start parameters object

                Returns:
                    WarmStartParams: copy of warm start parameters
            )pbdoc")
    ;

    py::class_<Interval>(m, "Interval",
            R"pbdoc(
                Interval class

                Implements interface from IntervalBase. This class owns its lower and upper bounds.
            )pbdoc")
        .def_readwrite("lb", &Interval::lb, "lower bound")
        .def_readwrite("ub", &Interval::ub, "upper bound")
        .def(py::init<zono_float, zono_float>(), py::arg("y_min"), py::arg("y_max"),
            R"pbdoc(
                Interval constructor.

                Args:
                    y_min (float): lower bound
                    y_max (float): upper bound
            )pbdoc")
        .def("__repr__", &Interval::print,
            R"pbdoc(
                print method for Interval

                Returns:
                    str: string representation of interval
            )pbdoc")
        .def("copy", &Interval::clone,
            R"pbdoc(
                Copy interval object

                Returns:
                    Interval: copy of interval
            )pbdoc")
        .def("__add__", &Interval::operator+, py::arg("other"),
            R"pbdoc(
                Interval addition

                Args:
                    other (Interval): rhs interval

                Returns:
                    Interval: self + other
            )pbdoc")
        .def("__sub__", &Interval::operator-, py::arg("other"),
            R"pbdoc(
                Interval subtraction

                Args:
                    other (Interval): rhs interval

                Returns:
                    Interval: self - other
            )pbdoc")
        .def("__mul__", [](const Interval& self, const Interval& other) -> Interval { return self*other; } ,
            py::arg("other"),
            R"pbdoc(
                Interval multiplication

                Args:
                    other (Interval): rhs interval

                Returns:
                    Interval: self * other
            )pbdoc")
        .def("__mul__", [](const Interval& self, const zono_float alpha) -> Interval { return self*alpha; },
            py::arg("alpha"),
            R"pbdoc(
                Interval multiplication with scalar

                Args:
                    alpha (float): scalar multiplier

                Returns:
                    Interval: alpha * self
            )pbdoc")
        .def("__truediv__", &Interval::operator/, py::arg("other"),
            R"pbdoc(
                Interval division

                Args:
                    other (Interval): rhs interval

                Returns:
                    Interval: self / other
            )pbdoc")
        .def("inv", &Interval::inv,
            R"pbdoc(
                Interval inverse

                Returns:
                    Interval: inverse of self
            )pbdoc")
        .def("intersect", &Interval::intersect, py::arg("other"),
            R"pbdoc(
                Interval intersection

                Args:
                    other (Interval): rhs interval

                Returns:
                    Interval: intersection of self and other
            )pbdoc")
        .def("is_empty", &Interval::is_empty,
            R"pbdoc(
                Checks whether interval is empty

                Returns:
                    bool: flag indicating whether interval is empty
            )pbdoc")
        .def("contains", &Interval::contains, py::arg("y"),
            R"pbdoc(
                Checks whether interval contains a value

                Args:
                    y (float): scalar value

                Returns:
                    bool: flag indicating if interval contains y
            )pbdoc")
        .def("is_single_valued", &Interval::is_single_valued,
            R"pbdoc(
                Checks whether interval is single-valued (i.e., width is 0 within numerical tolerance)

                Returns:
                    bool: flag indicating if interval is single-value
            )pbdoc")
        .def("width", &Interval::width,
            R"pbdoc(
                Gets width of interval (ub - lb)

                Returns:
                    float: width of interval
            )pbdoc")
        .def("center", &Interval::center,
            R"pbdoc(
                Gets center of interval (ub + lb) / 2

                Returns:
                    float: center of interval
            )pbdoc")
        .def("radius", &Interval::radius,
            R"pbdoc(
                Gets radius of interval

                Returns interval centered at zero with width equal to the width of the original interval

                Returns:
                    Interval: radius of interval
            )pbdoc")
        .def("sin", &Interval::sin,
            R"pbdoc(
                Compute interval containing sin(x) for all x in interval

                Returns:
                    Interval: interval containing sin(x)
            )pbdoc")
        .def("cos", &Interval::cos,
            R"pbdoc(
                Compute interval containing cos(x) for all x in interval

                Returns:
                    Interval: interval containing cos(x)
            )pbdoc")
        .def("tan", &Interval::tan,
            R"pbdoc(
                Compute interval containing tan(x) for all x in interval

                Returns:
                    Interval: interval containing tan(x)
            )pbdoc")
        .def("arcsin", &Interval::arcsin,
            R"pbdoc(
                Compute interval containing arcsin(x) for all x in interval

                Returns:
                    Interval: interval containing arcsin(x)
            )pbdoc")
        .def("arccos", &Interval::arccos,
            R"pbdoc(
                Compute interval containing arccos(x) for all x in interval

                Returns:
                    Interval: interval containing arccos(x)
            )pbdoc")
        .def("arctan", &Interval::arctan,
            R"pbdoc(
                Compute interval containing arctan(x) for all x in interval

                Returns:
                    Interval: interval containing arctan(x)
            )pbdoc")
        .def("exp", &Interval::exp,
            R"pbdoc(
                Compute interval containing exp(x) for all x in interval

                Returns:
                    Interval: interval containing exp(x)
            )pbdoc")
    ;

    py::class_<Box>(m, "Box", "Box (i.e., interval vector) class")
        .def(py::init<const Eigen::Vector<zono_float, -1>&, const Eigen::Vector<zono_float, -1>&>(), py::arg("x_lb"), py::arg("x_ub"),
            R"pbdoc(
                Constructor from intervals of lower and upper bounds

                Args:
                    x_lb (numpy.array): vector of lower bounds
                    x_ub (numpy.array): vector of upper bounds
            )pbdoc")
        .def("__repr__", &Box::print,
            R"pbdoc(
                print method for Box

                Returns:
                    str: string representation of Box
            )pbdoc")
        .def("__setitem__", [](Box& self, const int i, const Interval& val) -> void
            { self[i] = val; }, py::arg("i"), py::arg("val"),
            R"pbdoc(
                Set indexed interval in box to specified value

                Args:
                    i (int): index
                    val (Interval): new interval for index i in Box
            )pbdoc")
        .def("__getitem__", [](const Box& self, const int i) -> Interval
            { return self[i]; }, py::arg("i"),
            R"pbdoc(
                Get interval at index i

                Args:
                    i (int): index

                Returns:
                    Interval: interval at index i in Box
            )pbdoc")
        .def("size", &Box::size,
            R"pbdoc(
                Get size of Box object

                Returns:
                    int: size of box
            )pbdoc")
        .def("project", &Box::project, py::arg("x"),
            R"pbdoc(
                Projects vector onto the Box (in place)

                Args:
                    x (numpy.array): vector to be projected
            )pbdoc")
        .def("copy", &Box::clone,
            R"pbdoc(
                Copies Box object

                Returns:
                    Box: copy of object
            )pbdoc")
        .def("lower", &Box::lower,
            R"pbdoc(
                Get reference to lower bounds

                Returns:
                    numpy.array: lower bounds
            )pbdoc")
        .def("upper", &Box::upper,
            R"pbdoc(
                Get reference to upper bounds

                Returns:
                    numpy.array: upper bounds
            )pbdoc")
        .def("width", &Box::width,
            R"pbdoc(
                Get width of box.

                Specifically, this returns the max width for any interval in the box

                Returns:
                    float: width of box
            )pbdoc")
        .def("radius", &Box::radius,
            R"pbdoc(
                Get radius of box

                Returns box with intervals centered at zero with width equal to the width of the original box

                Returns:
                    Box: radius of box
            )pbdoc")
        .def("center", &Box::center,
            R"pbdoc(
                Gets center of box (x_ub + x_lb) / 2

                Returns:
                    numpy.array: center of interval
            )pbdoc")
        .def("contract", &Box::contract, py::arg("A"), py::arg("b"), py::arg("iter"),
            R"pbdoc(
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
            )pbdoc")
        .def("linear_map", [](const Box& self, const Eigen::SparseMatrix<zono_float, Eigen::RowMajor>& A) -> Box
            { return self.linear_map(A); }, py::arg("A"),
            R"pbdoc(
                Linear map of box based on interval arithmetic

                Args:
                    A (scipy.sparse.csr_matrix): linear map matrix

                Returns:
                    Box: linear mapped box
            )pbdoc")
        .def("dot", &Box::dot, py::arg("x"),
            R"pbdoc(
                Linear map with vector

                Args:
                    x (numpy.array): vector

                Returns:
                    Interval: result of linear map of box with vector
            )pbdoc")
        .def("__add__", &Box::operator+, py::arg("other"),
            R"pbdoc(
                Elementwise addition

                Args:
                    other (Box): rhs box

                Returns:
                    Box: self + other (elementwise)
            )pbdoc")
        .def("__sub__", &Box::operator-, py::arg("other"),
            R"pbdoc(
                Elementwise subtraction

                Args:
                    other (Box): rhs box

                Returns:
                    Box: self - other (elementwise)
            )pbdoc")
        .def("__mul__", [](const Box& self, const Box& other) -> Box { return self*other; } ,
            py::arg("other"),
            R"pbdoc(
                Elementwise multiplication

                Args:
                    other (Box): rhs box

                Returns:
                    Box: self * other (elementwise)
            )pbdoc")
        .def("__mul__", [](const Box& self, const zono_float alpha) -> Box { return self*alpha; },
            py::arg("alpha"),
            R"pbdoc(
                Elementwise multiplication with scalar

                Args:
                    alpha (float): scalar multiplier

                Returns:
                    Box: alpha * self (elementwise)
            )pbdoc")
        .def("__truediv__", &Box::operator/, py::arg("other"),
            R"pbdoc(
                Elementwise division

                Args:
                    other (Box): rhs box

                Returns:
                    Box: self / other (elementwise)
            )pbdoc")
    ;

    // interval matrix
    py::class_<IntervalMatrix>(m, "IntervalMatrix", "Interval matrix class")
        .def(py::init<const Eigen::SparseMatrix<zono_float>&, const Eigen::SparseMatrix<zono_float>&>(),
            py::arg("mat_lb"), py::arg("mat_ub"),
            R"pbdoc(
                IntervalMatrix constructor

                Args:
                    mat_lb (scipy.sparse.csc_matrix): matrix of lower bounds
                    mat_ub (scipy.sparse.csc_matrix): matrix of upper bounds
            )pbdoc")
        .def("center", &IntervalMatrix::center,
            R"pbdoc(
                Get center matrix

                Each element of center matrix is the center of the corresponding interval in the interval matrix

                Returns:
                    scipy.sparse.csc_matrix: center matrix
            )pbdoc")
        .def("diam", &IntervalMatrix::diam,
            R"pbdoc(
                Get diameter matrix

                Each element of the diameter matrix is the width of the corresponding interval in the interval matrix

                Returns:
                    scipy.sparse.csc_matrix: diameter matrix
            )pbdoc")
        .def("width", &IntervalMatrix::width,
            R"pbdoc(
                Get width of interval matrix

                Specifically, this returns the max width for any interval in the interval matrix

                Returns:
                    float: width of interval matrix
            )pbdoc")
        .def("radius", &IntervalMatrix::radius,
            R"pbdoc(
                Get radius of interval matrix

                Returns interval matrix with intervals centered at zero with width equal to the width of the original interval matrix

                Returns:
                    IntervalMatrix: radius of interval matrix
            )pbdoc")
        .def("__mul__", [](const IntervalMatrix& self, const Eigen::Vector<zono_float, -1>& v) -> Box
            { return self*v; }, py::arg("v"),
            R"pbdoc(
                IntervalMatrix multiplication with vector

                Args:
                    v (numpy.array): rhs vector

                Returns:
                    Box: resulting box
            )pbdoc")
        .def("__mul__", [](const IntervalMatrix& self, const Box& box) -> Box
            { return self*box; }, py::arg("box"),
            R"pbdoc(
                IntervalMatrix multiplication with Box

                Args:
                    box (Box): rhs box

                Returns:
                    Box: resulting box
            )pbdoc")
        .def("__mul__", [](const IntervalMatrix& self, const Eigen::SparseMatrix<zono_float, Eigen::RowMajor>& A) -> IntervalMatrix
            { return self*A; }, py::arg("A"),
            R"pbdoc(
                IntervalMatrix multiplication with matrix

                Args:
                    A (scipy.sparse.csr_matrix): rhs matrix

                Returns:
                    IntervalMatrix: resulting interval matrix
            )pbdoc")
        .def("__mul__", [](const IntervalMatrix& self, const IntervalMatrix& other) -> IntervalMatrix
            { return self*other; }, py::arg("other"),
            R"pbdoc(
                IntervalMatrix multiplication with another IntervalMatrix

                Args:
                    other (IntervalMatrix): rhs interval matrix

                Returns:
                    IntervalMatrix: resulting interval matrix
            )pbdoc")
        .def("__add__", &IntervalMatrix::operator+, py::arg("other"),
            R"pbdoc(
                IntervalMatrix addition

                Args:
                    other (IntervalMatrix): rhs interval matrix

                Returns:
                    IntervalMatrix: resulting interval matrix
            )pbdoc")
        .def("__sub__", &IntervalMatrix::operator-, py::arg("other"),
            R"pbdoc(
                IntervalMatrix subtraction

                Args:
                    other (IntervalMatrix): rhs interval matrix

                Returns:
                    IntervalMatrix: resulting interval matrix
            )pbdoc")
        .def("rows", &IntervalMatrix::rows,
            R"pbdoc(
                Get number of rows

                Returns:
                    int: number of rows
            )pbdoc")
        .def("cols", &IntervalMatrix::cols,
            R"pbdoc(
                Get number of columns

                Returns:
                    int: number of cols
            )pbdoc")
        .def("__repr__", &IntervalMatrix::print,
            R"pbdoc(
                Print method

                Returns:
                    str: string display of IntervalMatrix
            )pbdoc")
    ;


    // hybzono class
    py::class_<HybZono, py::smart_holder>(m, "HybZono", R"pbdoc(
            Hybrid zonotope class
             
            A hybrid zonotope is defined as:
            Z = {Gc * xi_c + Gb * xi_b + c | Ac * xi_c + Ab * xi_b = b, xi_c in [-1, 1]^nGc, xi_b in {-1, 1}^nGb}.
            Equivalently, the following shorthand can be used: Z = <Gc, Gb, c, Ac, Ab, b>.
            Optionally, in 0-1 form, the factors are xi_c in [0, 1]^nGc, xi_b in {0, 1}^nGb. 
            The set dimension is n, and the number of equality constraints is nC.
        )pbdoc")
        .def(py::init<const Eigen::SparseMatrix<zono_float>&, const Eigen::SparseMatrix<zono_float>&, const Eigen::Vector<zono_float, -1>&,
            const Eigen::SparseMatrix<zono_float>&, const Eigen::SparseMatrix<zono_float>&, const Eigen::Vector<zono_float, -1>&,
            bool, bool>(), py::arg("Gc"), py::arg("Gb"), py::arg("c"), py::arg("Ac"), py::arg("Ab"),
            py::arg("b"), py::arg("zero_one_form")=false, py::arg("sharp")=false, 
            R"pbdoc(
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
            )pbdoc")
        .def("set", &HybZono::set, py::arg("Gc"), py::arg("Gb"), py::arg("c"), py::arg("Ac"), py::arg("Ab"), 
            py::arg("b"), py::arg("zero_one_form")=false, py::arg("sharp")=false, 
            R"pbdoc(
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
            )pbdoc")
        .def("get_n", &HybZono::get_n,
            R"pbdoc(
                Returns dimension of set
                
                Returns:
                    int: n)
            )pbdoc")
        .def("get_nC", &HybZono::get_nC,
            R"pbdoc(
                Returns number of constraints in set definition
                
                Returns:
                    int: nC
            )pbdoc")
        .def("get_nG", &HybZono::get_nG,
            R"pbdoc(
                Returns number of generators in set definition
                
                Returns:
                    int: nG
            )pbdoc")
        .def("is_0_1_form", &HybZono::is_0_1_form, 
            R"pbdoc(
                Returns true if factors are in range [0,1], false if they are in range [-1,1].
                
                Returns:
                    bool: zero_one_form flag
            )pbdoc")
        .def("is_sharp", &HybZono::is_sharp, 
            R"pbdoc(
                Returns true if set is known to be sharp
                
                Returns:
                    bool: sharp flag
            )pbdoc")
        .def("get_G", &HybZono::get_G, 
            R"pbdoc(
                Returns generator matrix
                
                Returns:
                    scipy.sparse.csc_matrix: G
            )pbdoc")
        .def("get_c", &HybZono::get_c,
            R"pbdoc(
                Returns center vector
                
                Returns:
                    numpy.array: c
            )pbdoc")
        .def("get_A", &HybZono::get_A,
            R"pbdoc(
                Returns constraint matrix
                
                Returns:
                    scipy.sparse.csc_matrix: A
            )pbdoc")
        .def("get_b", &HybZono::get_b,
            R"pbdoc(
                Returns constraint vector
                
                Returns:
                    numpy.array: b
            )pbdoc")
        .def("get_nGc", &HybZono::get_nGc,
            R"pbdoc(
                Returns number of continuous generators in set definition
                
                Returns:
                    int: nGc
            )pbdoc")
        .def("get_nGb", &HybZono::get_nGb,
            R"pbdoc(
                Returns number of binary generators in set definition
                
                Returns:
                    int: nGb
            )pbdoc")
        .def("get_Gc", &HybZono::get_Gc,
            R"pbdoc(
                Returns continuous generator matrix
                
                Returns:
                    scipy.sparse.csc_matrix: Gc
            )pbdoc")
        .def("get_Gb", &HybZono::get_Gb,
            R"pbdoc(
                Returns binary generator matrix
                
                Returns:
                    scipy.sparse.csc_matrix: Gb
            )pbdoc")
        .def("get_Ac", &HybZono::get_Ac,
            R"pbdoc(
                Returns continuous constraint matrix
                
                Returns:
                    scipy.sparse.csc_matrix: Ac
            )pbdoc")
        .def("get_Ab", &HybZono::get_Ab,
            R"pbdoc(
                Returns binary constraint matrix
                
                Returns:
                    scipy.sparse.csc_matrix: Ab
            )pbdoc")
        .def("convert_form", &HybZono::convert_form, 
            R"pbdoc(
                Converts the set representation between -1-1 and 0-1 forms.
                
                This method converts the set representation between -1-1 and 0-1 forms. 
                If the set is in -1-1 form, then xi_c in [-1,1] and xi_b in {-1,1}.
                If the set is in 0-1 form, then xi_c in [0,1] and xi_b in {0,1}.
            )pbdoc")
        .def("remove_redundancy", &HybZono::remove_redundancy, py::arg("contractor_iter")=10,
            R"pbdoc(
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
            )pbdoc")
        .def("is_point", &HybZono::is_point,
            R"pbdoc(
                Polymorphic type checking
                
                Returns:
                    bool: true if set is a point
            )pbdoc")
        .def("is_zono", &HybZono::is_zono,
            R"pbdoc(
                Polymorphic type checking
                
                Returns:
                    bool: true if set is a zonotope
            )pbdoc")
        .def("is_conzono", &HybZono::is_conzono,
            R"pbdoc(
                Polymorphic type checking
                
                Returns:
                    bool: true if set is a constrained zonotope
            )pbdoc")
        .def("is_hybzono", &HybZono::is_hybzono,
            R"pbdoc(
                Polymorphic type checking
                
                Returns:
                    bool: true if set is a hybrid zonotope
            )pbdoc")
        .def("is_empty_set", &HybZono::is_empty_set,
            R"pbdoc(
                Polymorphic type checking

                Returns:
                    bool: true if set is a empty set object
            )pbdoc")
        .def("__repr__", &HybZono::print, 
            R"pbdoc(
                Returns set information as a string
                
                Returns:
                    str: set information
            )pbdoc")
        .def("optimize_over", [](const HybZono& self, const Eigen::SparseMatrix<zono_float> &P,
            const Eigen::Vector<zono_float, -1> &q, zono_float c,
            const OptSettings &settings, OptSolution* solution,
            const WarmStartParams& warm_start_params)-> Eigen::Vector<zono_float, -1>
            {
                auto sol_shared = std::make_shared<OptSolution>();
                auto z = self.optimize_over(P, q, c, settings, &sol_shared, warm_start_params);
                if (solution)
                    *solution = *sol_shared;
                return z;
            },
            "optimize over", py::arg("P"), py::arg("q"), py::arg("c")=0,
            py::arg("settings")=OptSettings(), py::arg("solution")=nullptr,
            py::arg("warm_start_params")=WarmStartParams(),
            R"pbdoc(
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
            )pbdoc")
        .def("project_point", [](const HybZono& self, const Eigen::Vector<zono_float, -1> &x,
            const OptSettings &settings, OptSolution* solution,
            const WarmStartParams& warm_start_params) -> Eigen::Vector<zono_float, -1>
            {
                auto sol_shared = std::make_shared<OptSolution>();
                auto z = self.project_point(x, settings, &sol_shared, warm_start_params);
                if (solution)
                    *solution = *sol_shared;
                return z;
            },
            py::arg("x"), py::arg("settings")=OptSettings(),
            py::arg("solution")=nullptr, py::arg("warm_start_params")=WarmStartParams(),
            R"pbdoc(
                Returns the projection of the point x onto the set object.
                
                Args:
                    x (numpy.array): point to be projected
                    settings (OptSettings, optional): optimization settings structure
                    solution (OptSolution, optional): optimization solution structure pointer, populated with result
                    warm_start_params (WarmStartParams, optional): warm start parameters structure

                Returns:
                    numpy.array: point z in the current set
            )pbdoc")
        .def("is_empty", [](const HybZono& self, const OptSettings &settings,
            OptSolution* solution, const WarmStartParams& warm_start_params) -> bool
            {
                auto sol_shared = std::make_shared<OptSolution>();
                const bool empty = self.is_empty(settings, &sol_shared, warm_start_params);
                if (solution)
                    *solution = *sol_shared;
                return empty;
            },
            py::arg("settings")=OptSettings(),
            py::arg("solution")=nullptr, py::arg("warm_start_params")=WarmStartParams(),
            R"pbdoc(
                Returns true if the set is provably empty, false otherwise.
                
                Args:
                    settings (OptSettings, optional): optimization settings structure
                    solution (OptSolution, optional): optimization solution structure pointer, populated with result
                    warm_start_params (WarmStartParams, optional): warm start parameters structure

                Returns:
                    bool: flag indicating whether set is provably empty
            )pbdoc")
        .def("support", [](HybZono& self, const Eigen::Vector<zono_float, -1> &d,
            const OptSettings &settings, OptSolution* solution,
            const WarmStartParams& warm_start_params) -> zono_float
            {
                auto sol_shared = std::make_shared<OptSolution>();
                const zono_float s = self.support(d, settings, &sol_shared, warm_start_params);
                if (solution)
                    *solution = *sol_shared;
                return s;
            },
            py::arg("d"), py::arg("settings")=OptSettings(),
            py::arg("solution")=nullptr, py::arg("warm_start_params")=WarmStartParams(),
            R"pbdoc(
                Computes support function of the set in the direction d.
                
                Args:
                    d (numpy.array): vector defining direction for support function
                    settings (OptSettings, optional): optimization settings structure
                    solution (OptSolution, optional): optimization solution structure pointer, populated with result
                    warm_start_params (WarmStartParams, optional): warm start parameters structure

                Returns:
                    float: support value

                Solves max_{z in Z} <z, d> where <., .> is the inner product
            )pbdoc")
        .def("contains_point", [](const HybZono& self, const Eigen::Vector<zono_float, -1> &x,
            const OptSettings &settings, OptSolution* solution,
            const WarmStartParams& warm_start_params) -> bool
            {
                auto sol_shared = std::make_shared<OptSolution>();
                const bool contains = self.contains_point(x, settings, &sol_shared, warm_start_params);
                if (solution)
                    *solution = *sol_shared;
                return contains;
            },
            py::arg("x"), py::arg("settings")=OptSettings(),
            py::arg("solution")=nullptr, py::arg("warm_start_params")=WarmStartParams(),
            R"pbdoc(
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
            )pbdoc")
        .def("bounding_box", [](HybZono& self,
            const OptSettings &settings, OptSolution* solution,
            const WarmStartParams& warm_start_params) -> Box
            {
                auto sol_shared = std::make_shared<OptSolution>();
                const Box bb = self.bounding_box(settings, &sol_shared, warm_start_params);
                if (solution)
                    *solution = *sol_shared;
                return bb;
            },
            py::arg("settings")=OptSettings(),
            py::arg("solution")=nullptr, py::arg("warm_start_params")=WarmStartParams(),
            R"pbdoc(
                Computes a bounding box of the set object as a Box object.
                
                Args:
                    settings (OptSettings, optional): optimization settings structure
                    solution (OptSolution, optional): optimization solution structure pointer, populated with result
                    warm_start_params (WarmStartParams, optional): warm start parameters structure

                Returns:
                    Box: bounding box of the set

                In general, solves 2*n support optimizations where n is the set dimension to compute a bounding box.
            )pbdoc")
        .def("convex_relaxation", &HybZono::convex_relaxation, 
            R"pbdoc(
                Computes the convex relaxation of the hybrid zonotope.
                
                Returns:
                    ConZono: Constrained zonotope Z = <[Gc, Gb], c, [Ac, Ab,], b>

                This method returns the convex relaxation of the hybrid zonotope.
                If the set is sharp, the convex relaxation is the convex hull.
            )pbdoc")
        .def("get_leaves", [](const HybZono& self, bool remove_redundancy,
            const OptSettings &settings, OptSolution* solution,
            int n_leaves, int contractor_iter) -> std::vector<ConZono>
            {
                auto sol_shared = std::make_shared<OptSolution>();
                const auto leaves = self.get_leaves(remove_redundancy, settings, &sol_shared, n_leaves, contractor_iter);
                if (solution)
                    *solution = *sol_shared;
                return leaves;
            },
            py::arg("remove_redundancy")=true, py::arg("settings")=OptSettings(), py::arg("solution")=nullptr,
            py::arg("n_leaves")=std::numeric_limits<int>::max(), py::arg("contractor_iter")=100, 
            R"pbdoc(
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
            )pbdoc")
        .def("complement", [](HybZono& self, zono_float delta_m,
            bool remove_redundancy, const OptSettings &settings, OptSolution* solution,
            int n_leaves, int contractor_iter) -> std::unique_ptr<HybZono>
            {
                auto sol_shared = std::make_shared<OptSolution>();
                auto Z = self.complement(delta_m, remove_redundancy, settings, &sol_shared, n_leaves, contractor_iter);
                if (solution)
                    *solution = *sol_shared;
                return Z;
            },
            py::arg("delta_m")=100, py::arg("remove_redundancy")=true, py::arg("settings")=OptSettings(),
            py::arg("solution")=nullptr, py::arg("n_leaves")=std::numeric_limits<int>::max(), py::arg("contractor_iter")=100,
            R"pbdoc(
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
            X = {G \xi + c | A \xi = b, \xi \in [-1-delta_m, 1+delta+m]^{nG}}.
            )pbdoc")
        .def("copy", &HybZono::clone,
            R"pbdoc(
            Creates a copy of the hybrid zonotope object.

            Returns:
                HybZono: A copy of the hybrid zonotope object.
            )pbdoc")
    ;

    // conzono class
    py::class_<ConZono, HybZono /* parent type */, py::smart_holder>(m, "ConZono",
            R"pbdoc(
                Constrained zonotope class
                
                A constrained zonotope is defined as:
                Z = {G \xi + c | A \xi = b, \xi in [-1, 1]^nG}.
                Equivalently, the following shorthand can be used: Z = <G, c, A, b>.
                Optionally, in 0-1 form, the factors are xi in [0,1].
                The set dimension is n, and the number of equality constraints is nC.
            )pbdoc")
        .def(py::init<const Eigen::SparseMatrix<zono_float>&, const Eigen::Vector<zono_float, -1>&,
            const Eigen::SparseMatrix<zono_float>&, const Eigen::Vector<zono_float, -1>&, 
            bool>(), py::arg("G"), py::arg("c"), py::arg("A"), py::arg("b"), py::arg("zero_one_form")=false,
            R"pbdoc(
                ConZono constructor
                
                Args:
                    G (scipy.sparse.csc_matrix): generator matrix
                    c (numpy.array): center
                    A (scipy.sparse.csc_matrix): constraint matrix
                    b (numpy.array): constraint vector
                    zero_one_form (bool, optional): true if set is in 0-1 form
            )pbdoc")
        .def("set", &ConZono::set, py::arg("G"), py::arg("c"), py::arg("A"), py::arg("b"), 
            py::arg("zero_one_form")=false,
            R"pbdoc(
                Reset constrained zonotope object with the given parameters.
                
                Args:
                    G (scipy.sparse.csc_matrix): generator matrix
                    c (numpy.array): center
                    A (scipy.sparse.csc_matrix): constraint matrix
                    b (numpy.array): constraint vector
                    zero_one_form (bool, optional): true if set is in 0-1 form
            )pbdoc")
        .def("to_zono_approx", &ConZono::to_zono_approx,
           R"pbdoc(
                Compute outer approximation of constrained zonotope as zonotope using SVD

                Returns:
                    Zono: Zonotope over-approximation
            )pbdoc")
        .def("constraint_reduction", &ConZono::constraint_reduction,
            R"pbdoc(
                Execute constraint reduction algorithm from Scott et. al. 2016

                Removes one constraint and one generator from the constrained zonotope.
                The resulting set is an over-approximation of the original set.
            )pbdoc")
    ;

    // zono class
    py::class_<Zono, ConZono /* parent type */, py::smart_holder>(m, "Zono",
            R"pbdoc(
                Zonotope class
                
                A zonotope is defined as:
                Z = {G \xi + c | \xi in [-1, 1]^nG}.
                Equivalently, the following shorthand can be used: Z = <G, c>.
                Optionally, in 0-1 form, the factors are xi in [0,1].
                The set dimension is n, and the number of generators is nG.
            )pbdoc")
        .def(py::init<const Eigen::SparseMatrix<zono_float>&, const Eigen::Vector<zono_float, -1>&, bool>(), 
            py::arg("G"), py::arg("c"), py::arg("zero_one_form")=false,
            R"pbdoc(
                Zono constructor
                
                Args:
                    G (scipy.sparse.csc_matrix): generator matrix
                    c (numpy.array): center
                    zero_one_form (bool, optional): true if set is in 0-1 form
            )pbdoc")
        .def("set", &Zono::set, "set zonotope", py::arg("G"), py::arg("c"), py::arg("zero_one_form")=false,
            R"pbdoc(
                Reset zonotope object with the given parameters.
                
                Args:
                    G (scipy.sparse.csc_matrix): generator matrix
                    c (numpy.array): center
                    zero_one_form (bool, optional): true if set is in 0-1 form
            )pbdoc")
        .def("reduce_order", &Zono::reduce_order, py::arg("n_o"),
            R"pbdoc(
                Perform zonotope order reduction.

                Args:
                    n_o (int): desired order, must be greater than or equal to the dimension of the set

                Returns:
                    zonotope with order n_o
            )pbdoc")
        .def("get_volume", &Zono::get_volume,
            R"pbdoc(
                    Get volume of zonotope.

                    Reference: Gover and Krikorian 2010, "Determinants and the volumes of parallelotopes and zonotopes"
                    Requires nG choose n determinant computations.

                    Returns:
                        float: volume of zonotope
            )pbdoc")
        .def("get_center", &Zono::get_center,
            R"pbdoc(
                Get center of zonotope.

                Returns:
                    numpy.array: center vector
            )pbdoc")
    ;

    // point class
    py::class_<Point, Zono /* parent type */, py::smart_holder>(m, "Point",
            R"pbdoc(
                Point class
                
                A point is defined entirely by the center vector c.
            )pbdoc")
        .def(py::init<const Eigen::Vector<zono_float, -1>&>(), py::arg("c"),
            R"pbdoc(
                Point constructor
                
                Args:
                    c (numpy.array): center vector
            )pbdoc")
        .def("set", &Point::set, "set point", py::arg("c"),
            R"pbdoc(
                Reset point object with the given parameters.
                
                Args:
                    c (numpy.array): center vector
            )pbdoc")
    ;

    // empty set class
    py::class_<EmptySet, ConZono, py::smart_holder>(m, "EmptySet",
            R"pbdoc(
                Empty Set class

                Used to facilitate set operations with trivial solutions when one of the sets is an empty set.
            )pbdoc")
        .def(py::init<int>(), py::arg("n"),
            R"pbdoc(
                EmptySet constructor

                Args:
                    n (int): dimension
            )pbdoc")
    ;

    // inequalities
    py::class_<IneqTerm>(m, "IneqTerm", "Structure containing term in 0-1 inequality.")
        .def(py::init<int, zono_float>(), "IneqTerm constructor", py::arg("idx"), py::arg("coeff"))
        .def_readwrite("idx", &IneqTerm::idx, "index of variable")
        .def_readwrite("coeff", &IneqTerm::coeff, "coefficient of variable")
        .def("copy", [](const IneqTerm &self) -> IneqTerm { return self; },
            "Creates a copy of the IneqTerm object.")
    ;

    py::enum_<IneqType>(m, "IneqType", "Enumeration to select inequality direction / use equality.")
        .value("LESS", LESS, "Strictly less than")
        .value("LESS_OR_EQUAL", LESS_OR_EQUAL, "Less than or equal to")
        .value("EQUAL", EQUAL, "Equal to")
        .value("GREATER_OR_EQUAL", GREATER_OR_EQUAL, "Greater than or equal to")
        .value("GREATER", GREATER, "Strictly greater than")
        .export_values()
    ;

    py::class_<Inequality>(m, "Inequality", "Inequality class")
        .def(py::init<int>(), py::arg("n_dims"),
            R"pbdoc(
                Constructs inequality, must specify number of dimensions.
                
                Args:
                    n_dims (int): number of dimensions in the inequality
            )pbdoc")
        .def("add_term", &Inequality::add_term, py::arg("idx"), py::arg("coeff"),
            R"pbdoc(
                Adds a term to the inequality.
                
                Args:
                    idx (int): index of variable in the inequality
                    coeff (float): coefficient of the variable
            )pbdoc")
        .def("set_rhs", &Inequality::set_rhs, py::arg("rhs"), 
            R"pbdoc(
                Sets right-hand side of the inequality.
                
                Args:
                    rhs (float): right-hand side value
            )pbdoc")
        .def("set_ineq_type", &Inequality::set_ineq_type, py::arg("type"), 
            R"pbdoc(
                Sets the direction of the inequality or sets it to be an equality
                
                Args:
                    type (IneqType): inequality type (e.g., less than or equal, greater than or equal, or equal)
            )pbdoc")
        .def("get_terms", &Inequality::get_terms, 
            R"pbdoc(
                Get reference to terms in the inequality.
                
                Returns:
                    list[IneqTerm]: reference to terms member
            )pbdoc")
        .def("get_rhs", &Inequality::get_rhs, 
            R"pbdoc(
                Get right-hand side of the inequality.
                
                Returns:
                    float: right-hand side value (rhs member)
            )pbdoc")
        .def("get_ineq_type", &Inequality::get_ineq_type, 
            R"pbdoc(
                Get inequality type / direction.
                
                Returns:
                    IneqType: inequality type (type member)
            )pbdoc")
        .def("get_n_dims", &Inequality::get_n_dims, 
            R"pbdoc(
                Get number of dimensions for inequality.
                
                Returns:
                    int: number of dimensions (n_dims member)
            )pbdoc")
        .def("copy", [](const Inequality &self) -> Inequality { return self; },
            R"pbdoc(
                Creates a copy of the Inequality object.

                Returns:
                    Inequality: A copy of the Inequality object.
            )pbdoc")
    ;

    // set operations
    m.def("affine_map", &affine_map, py::arg("Z"), py::arg("R"), py::arg("s")=Eigen::Vector<zono_float, -1>(),
        R"pbdoc(
            Returns affine map R*Z + s of set Z
            
            Args:
                Z (HybZono): zonotopic set
                R (scipy.sparse.csc_matrix): affine map matrix
                s (numpy.array, optional): vector offset
            
            Returns:
                HybZono: zonotopic set
        )pbdoc");
    m.def("affine_inclusion", &affine_inclusion, py::arg("Z"), py::arg("R"), py::arg("s")=Eigen::Vector<zono_float, -1>(),
        R"pbdoc(
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
        )pbdoc");
    m.def("project_onto_dims", &project_onto_dims, py::arg("Z"), py::arg("dims"),
        R"pbdoc(
            Projects set Z onto the dimensions specified in dims.
            
            Args:
                Z (HybZono): zonotopic set
                dims (list[int]): list of dimensions to project onto
            
            Returns:
                HybZono: zonotopic set
        )pbdoc");
    m.def("minkowski_sum", &minkowski_sum, py::arg("Z1"), py::arg("Z2"),
        R"pbdoc(
            Computes Minkowski sum of two sets Z1 and Z2.
            
            Args:
                Z1 (HybZono): zonotopic set
                Z2 (HybZono): zonotopic set
            
            Returns:
                HybZono: zonotopic set
        )pbdoc");
    m.def("pontry_diff", &pontry_diff, py::arg("Z1"), py::arg("Z2"), py::arg("exact")=false,
        R"pbdoc(
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
        )pbdoc");
    m.def("intersection", &intersection, py::arg("Z1"), py::arg("Z2"), py::arg("R")=Eigen::SparseMatrix<zono_float>(),
        R"pbdoc(
            Computes the generalized intersection of sets Z1 and Z2 over the matrix R.
            
            Args:
                Z1 (HybZono): zonotopic set
                Z2 (HybZono): zonotopic set
                R (scipy.sparse.csc_matrix, optional): affine map matrix
            
            Returns:
                HybZono: zonotopic set
        )pbdoc");
    m.def("intersection_over_dims", &intersection_over_dims, py::arg("Z1"), py::arg("Z2"), py::arg("dims"),
        R"pbdoc(
            Computes the intersection of sets Z1 and Z2 over the specified dimensions.
            
            Args:
                Z1 (HybZono): zonotopic set
                Z2 (HybZono): zonotopic set
                dims (list[int]): list of dimensions
            
            Returns:
                HybZono: zonotopic set
        )pbdoc");
    m.def("halfspace_intersection", &halfspace_intersection, py::arg("Z"), py::arg("H"), py::arg("f"), py::arg("R")=Eigen::SparseMatrix<zono_float>(),
        R"pbdoc(
            Computes the intersection generalized intersection of set Z with halfspace H*x <= f over matrix R.
            
            Args:
                Z (HybZono): zonotopic set
                H (scipy.sparse.csc_matrix): halfspace matrix
                f (numpy.array): halfspace vector
                R (scipy.sparse.csc_matrix, optional): affine map matrix
            
            Returns:
                HybZono: zonotopic set
        )pbdoc");
    m.def("union_of_many", &union_of_many, py::arg("Z_list"), py::arg("preserve_sharpness")=false, py::arg("expose_indicators")=false,
        R"pbdoc(
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
        )pbdoc");
    m.def("convex_hull", &convex_hull, py::arg("Zs"),
        R"pbdoc(
            Computes the convex hull of several sets

            Computes convex hull of sets {Z0, Z1, ..., Zn}.
            If Zi is a hybrid zonotope, it must be sharp or this function will throw an error.

            Args:
                Zs (list[HybZono]): sets for which convex hull is to be computed.

            Returns:
                ConZono: constrained zonotop convex hull
        )pbdoc");
    m.def("cartesian_product", &cartesian_product, py::arg("Z1"), py::arg("Z2"),
        R"pbdoc(
            Computes the Cartesian product of two sets Z1 and Z2.
            
            Args:
                Z1 (HybZono): zonotopic set
                Z2 (HybZono): zonotopic set
            
            Returns:
                HybZono: zonotopic set
        )pbdoc");
    m.def("constrain", &constrain, py::arg("Z"), py::arg("ineqs"), py::arg("R")=Eigen::SparseMatrix<zono_float>(),
        R"pbdoc(
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
        )pbdoc");
    m.def("set_diff", [](const HybZono& Z1, HybZono& Z2, zono_float delta_m,
            bool remove_redundancy, const OptSettings &settings, OptSolution* solution,
            int n_leaves, int contractor_iter) -> std::unique_ptr<HybZono>
            {
                auto sol_shared = std::make_shared<OptSolution>();
                auto Z = set_diff(Z1, Z2, delta_m, remove_redundancy, settings, &sol_shared, n_leaves, contractor_iter);
                if (solution)
                    *solution = *sol_shared;
                return Z;
            },
            py::arg("Z1"), py::arg("Z2"), py::arg("delta_m")=100, py::arg("remove_redundancy")=true,
            py::arg("settings")=OptSettings(), py::arg("solution")=nullptr, py::arg("n_leaves")=std::numeric_limits<int>::max(), py::arg("contractor_iter")=10,
            R"pbdoc(
            Set difference Z1 \\ Z2

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
        )pbdoc");

    // global setup functions
    m.def("vrep_2_conzono", &vrep_2_conzono, py::arg("Vpoly"), 
        R"pbdoc(
            Builds a constrained zonotope from a vertex representation polytope.

            Args:
                Vpoly (numpy.array): vertices of V-rep polytope
            
            Returns:
                ConZono: constrained zonotope
            
            Vpoly is a matrix where each row is a vertex of the polytope.
        )pbdoc");
    m.def("vrep_2_hybzono", &vrep_2_hybzono, py::arg("Vpolys"), py::arg("expose_indicators")=false,
        R"pbdoc(
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
        )pbdoc");
    m.def("zono_union_2_hybzono", &zono_union_2_hybzono, py::arg("Zs"), py::arg("expose_indicators")=false,
        R"pbdoc(
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
        )pbdoc");
    m.def("make_regular_zono_2D", &make_regular_zono_2D, py::arg("radius"), py::arg("n_sides"), py::arg("outer_approx")=false, py::arg("c")=Eigen::Vector<zono_float, 2>::Zero(),
        R"pbdoc(
            Builds a 2D regular zonotope with a given radius and number of sides.
            
            Args:
                radius (float): radius of the zonotope
                n_sides (int): number of sides (must be an even number >= 4)
                outer_approx (bool, optional): flag to do an outer approximation instead of an inner approximation
                c (numpy.array, optional): center vector
            
            Returns:
                Zono: zonotope
        )pbdoc");
    m.def("interval_2_zono", &interval_2_zono, py::arg("box"),
        R"pbdoc(
            Builds a zonotope from a Box object.
            
            Args:
                box (Box): Box object (vector of intervals)
            
            Returns:
                Zono: zonotope
        )pbdoc");
}