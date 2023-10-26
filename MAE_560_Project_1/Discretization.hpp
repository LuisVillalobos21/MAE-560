#pragma once
#include <Eigen/Dense>

struct Discretization
{
    Eigen::MatrixXd laplacian_matrix;

    Discretization(const int number_points);

    Eigen::MatrixXd createLaplacianOperator(const int number_points);
};