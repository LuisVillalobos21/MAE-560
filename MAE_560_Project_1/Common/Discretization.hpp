#pragma once
#include <Eigen/Dense>
#include "Parameters.hpp"
#include "Mesh.hpp"

struct Discretization
{
    Eigen::MatrixXd laplacian_matrix;

    Discretization(const int number_points);

    Eigen::MatrixXd createLaplacianOperator(const int number_points);
};

struct CNDiscretization
{
    Eigen::MatrixXd laplacian_matrix;
    Eigen::MatrixXd cn_matrix;
    CNDiscretization(const int number_points, const Parameters& params, const Mesh& mesh);
    Eigen::MatrixXd createLaplacianOperator(const int number_points);
    Eigen::MatrixXd createCNMatrix(const Parameters& params, const Mesh& mesh);
};