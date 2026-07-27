#include "Eigen/Dense"
#include "Eigen/Sparse"

#include "zonoopt/CholeskyUtilities.hpp"
#include "zonoopt/Defines.hpp"

namespace ZonoOpt::detail
{
    void get_LDLT_data(const Eigen::SimplicialLDLT<Eigen::SparseMatrix<zono_float>>& solver, LDLT_data& data)
    {
        data.L = solver.matrixL();
        data.Dinv = solver.vectorD().cwiseInverse().asDiagonal();
        data.P = solver.permutationP();
        data.Pinv = solver.permutationPinv();
        data.factorized = true;
    }

    Eigen::Vector<zono_float, -1> solve_LDLT(const LDLT_data& data, const Eigen::Vector<zono_float, -1>& b)
    {
        const Eigen::Vector<zono_float, -1> bbar = data.P * b;
        const Eigen::Vector<zono_float, -1> y = data.Dinv * data.L.template triangularView<Eigen::Lower>().solve(bbar);
        const Eigen::Vector<zono_float, -1> xbar = data.L.transpose().template triangularView<Eigen::Upper>().solve(y);
        return data.Pinv * xbar;
    }

    void affine_set_projection(Eigen::Ref<Eigen::Vector<zono_float, -1>> z, const Eigen::SparseMatrix<zono_float>& A,
                               const Eigen::SparseMatrix<zono_float>& AT, const Eigen::Vector<zono_float, -1>& b,
                               const LDLT_data& AAT_ldlt)
    {
        const Eigen::Vector<zono_float, -1> bmAz = b - A * z;
        const Eigen::Vector<zono_float, -1> bmAz_sol = solve_LDLT(AAT_ldlt, bmAz);
        z += AT * bmAz_sol;
    }
} // end namespace detail
// end namespace ZonoOpt
