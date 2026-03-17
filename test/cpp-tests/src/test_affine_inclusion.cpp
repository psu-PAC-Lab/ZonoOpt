#include "ZonoOpt.hpp"
#include "unit_test_utilities.hpp"
#include <cstdlib>
#include <iostream>

using namespace ZonoOpt;

std::tuple<Zono, Eigen::MatrixXd, Eigen::MatrixXd, Eigen::VectorXd> random_example(int n, int nG, int n_out)
{
    // generate zonotope
    Eigen::MatrixXd G = Eigen::MatrixXd::Random(n, nG);
    Eigen::VectorXd c = Eigen::VectorXd::Random(n);
    auto Z = Zono(G.sparseView(), c);

    // offset
    Eigen::VectorXd s = 10.*Eigen::VectorXd::Random(n_out);

    // random affine map
    Eigen::MatrixXd R = Eigen::MatrixXd::Random(n_out, n);
    auto R_lb = R;
    auto R_ub = R;
    for (int i=0; i<n_out; ++i)
    {
        for (int j=0; j<n; ++j)
        {
            double drij = 0.1*std::rand()/RAND_MAX;
            R_lb(i,j) = R(i,j) - drij;
            R_ub(i,j) = R(i,j) + drij;
        }
    }

    return std::make_tuple(Z, R_lb, R_ub, s);
}


int main()
{
    srand(0);

    for (int seed=0; seed<10; ++seed)
    {
        auto [Z, R_lb, R_ub, s] = random_example(5, 10, 3);

        // affine inclusion
        IntervalMatrix Rint (R_lb, R_ub);
        auto Z_inc = affine_inclusion(Z, Rint, s);

        // lower and upper bounds
        auto Z_lb = affine_map(Z, R_lb.sparseView(), s);
        auto Z_ub = affine_map(Z, R_ub.sparseView(), s);

        // make sure centers of Z_lb and Z_ub are inside Z_inc
        Zono Z_lb_zono = *dynamic_cast<Zono*>(Z_lb.get());
        Zono Z_ub_zono = *dynamic_cast<Zono*>(Z_ub.get());

        std::stringstream err_ss;
        err_ss << "Zonotope " << *Z_inc << " does not contain point " << Z_lb_zono.get_center();

        test_assert(Z_inc->contains_point(Z_lb_zono.get_center()), err_ss.str());

        err_ss.str("");
        err_ss << "Zonotope " << *Z_inc << " does not contain point " << Z_ub_zono.get_center();

        test_assert(Z_inc->contains_point(Z_ub_zono.get_center()), err_ss.str());
    }

    return 0;
}