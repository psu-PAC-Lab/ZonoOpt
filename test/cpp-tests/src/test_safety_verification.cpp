#include "ZonoOpt.hpp"
#include "unit_test_utilities.hpp"

int main()
{
    // System dynamics
    constexpr double dt = 0.1;
    Eigen::Matrix<double, 2, 2> A;
    A << 1, dt,
         0, 1;
    Eigen::Matrix<double, 2, 1> B;
    B << 0.5*dt*dt,
         dt;

    // Initial set: box [-1.0, 1.0] x [-0.1, 0.1]
    Eigen::Vector2d x0_min, x0_max;
    x0_min << -1.0, -0.1;
    x0_max <<  1.0,  0.1;
    ZonoOpt::ZonoPtr X0 = ZonoOpt::interval_2_zono(ZonoOpt::Box(x0_min, x0_max));

    // Input set: box [-0.2, 0.2]
    Eigen::Vector<double, 1> u_min, u_max;
    u_min << -0.2;
    u_max <<  0.2;
    ZonoOpt::ZonoPtr U = ZonoOpt::interval_2_zono(ZonoOpt::Box(u_min, u_max));

    // Disturbance set: affine map of octagon
    ZonoOpt::ZonoPtr W = ZonoOpt::make_regular_zono_2D(1.0, 8);
    Eigen::SparseMatrix<double> W_map(2, 2);
    W_map.insert(0, 0) = 0.01;
    W_map.insert(1, 1) = 0.05;
    W = ZonoOpt::affine_map(*W, W_map);

    // Compute reachable set over 10 time steps
    ZonoOpt::ZonoPtr X = std::move(X0);
    for (int k=0; k<10; ++k)
    {
        X = ZonoOpt::affine_map(*X, A.sparseView());
        X = ZonoOpt::minkowski_sum(*X, *ZonoOpt::affine_map(*U, B.sparseView()));
        X = ZonoOpt::minkowski_sum(*X, *W);
    }

    // Unsafe set
    Eigen::Matrix<double, 4, 2> verts;
    verts << 1.3, 0.0,
             1.6, 0.8,
             2.0, -0.4,
             2.3, 0.6;
    const ZonoOpt::ZonoPtr O = ZonoOpt::vrep_2_conzono(verts);

    // expect that intersection of X and O is empty
    test_assert(ZonoOpt::intersection(*X, *O)->is_empty(), "Expected intersection of X and O to be empty");

    return 0;
}