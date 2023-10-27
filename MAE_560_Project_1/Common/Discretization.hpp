#pragma once
#include <Eigen/Dense>
#include "Parameters.hpp"
#include "Mesh.hpp"

Eigen::MatrixXd createLaplacianMatrix(const int number_points);

Eigen::MatrixXd createDifferentiationMatrix(const int number_points);

struct Discretization
{
    Eigen::MatrixXd laplacian_matrix;

    Discretization(const int number_points);
};

struct CNDiscretization
{
    Eigen::MatrixXd laplacian_matrix;
    Eigen::MatrixXd cn_matrix;

    CNDiscretization(const int number_points, const Parameters& params, const Mesh& mesh);
    Eigen::MatrixXd createCNMatrix(const Parameters& params, const Mesh& mesh);
};

struct LWDiscretization
{
    Eigen::MatrixXd laplacian_matrix;
    Eigen::MatrixXd differentiaion_matrix;

    LWDiscretization(const int number_points);
};