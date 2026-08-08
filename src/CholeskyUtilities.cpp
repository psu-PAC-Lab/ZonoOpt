#include "Eigen/Dense"
#include "Eigen/Sparse"

#include "zonoopt/CholeskyUtilities.hpp"
#include "zonoopt/Defines.hpp"

namespace ZonoOpt::detail
{
    void affine_set_projection(Eigen::Ref<Eigen::Vector<zono_float, -1>> z, const Eigen::SparseMatrix<zono_float>& A,
                               const Eigen::SparseMatrix<zono_float>& AT, const Eigen::Vector<zono_float, -1>& b,
                               const LDLT_solver& AAT_ldlt)
    {
        const Eigen::Vector<zono_float, -1> bmAz = b - A * z;
        const Eigen::Vector<zono_float, -1> bmAz_sol = AAT_ldlt.solve(bmAz);
        z += AT * bmAz_sol;
    }
} // end namespace detail
// end namespace ZonoOpt
