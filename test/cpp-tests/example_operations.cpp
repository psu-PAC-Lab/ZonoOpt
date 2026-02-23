#include <iostream>
#include <exception>

// #define zono_double double
#include "ZonoOpt.hpp"
using namespace ZonoOpt;

int main()
{
    // settings struct
    OptSettings settings;
    settings.verbose = true;
    settings.n_threads_bnb = 1;

    // create conzono
    Eigen::SparseMatrix<double> G(2, 2);
    G.insert(0, 0) = 1;
    G.insert(1, 1) = 1;

    Eigen::Vector<double, -1> c(2);
    c.setZero();

    Eigen::SparseMatrix<double> A(1, 2);
    A.insert(0, 0) = 1;

    Eigen::Vector<double, -1> b(1);
    b(0) = 1;

    ConZono Z1 (G, c, A, b, true);
    Zono Z2 (G, c, false);

    // try minkowski sum
    ZonoPtr Z_sum = minkowski_sum(Z1, Z2);
    std::cout << "Z_sum is conzono?: " << Z_sum->is_conzono() << std::endl;
    std::cout << "Z_sum: " << *Z_sum << std::endl;

    // union
    ZonoPtr U = union_of_many({std::shared_ptr<HybZono>(Z1.clone()), std::shared_ptr<HybZono>(Z2.clone())});
    std::cout << "U: " << *U << std::endl;
    std::cout << "U is empty? " << U->is_empty(settings) << std::endl;

    // convex relaxation
    auto C = U->convex_relaxation();
    std::cout << "C: " << *C << std::endl;

    // project point
    Eigen::VectorXd x (C->get_n());
    x.setOnes();

    std::cout << "Point projection onto C: " << std::endl;
    std::cout << C->project_point(x, settings) << std::endl;

    std::cout << "Point projection onto U: " << std::endl;
    std::cout << U->project_point(x, settings) << std::endl;


    // build hybzono with redundant constraints
    Eigen::SparseMatrix<double> Gc_h(2, 2);
    Gc_h.setIdentity();
    Eigen::SparseMatrix<double> Gb_h(2, 2);
    Gb_h.setIdentity();
    Eigen::Vector<double, 2> c_h;
    c_h.setZero();

    Eigen::MatrixXd Acd_h = Eigen::MatrixXd::Ones(2, 2);
    Eigen::MatrixXd Abd_h = Eigen::MatrixXd::Ones(2, 2);
    Eigen::SparseMatrix<double> Ac_h = Acd_h.sparseView();
    Eigen::SparseMatrix<double> Ab_h = Abd_h.sparseView();
    Eigen::Vector<double, 2> b_h;
    b_h << 1, 1;

    HybZono Zh (Gc_h, Gb_h, c_h, Ac_h, Ab_h, b_h, true);
    Zh.remove_redundancy();
    
    std::cout << "Zh: " << Zh << std::endl;

    return 0;
}