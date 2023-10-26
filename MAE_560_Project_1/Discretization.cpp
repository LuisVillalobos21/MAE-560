#include "Discretization.hpp"

Discretization::Discretization(const int number_points)
{
    laplacian_matrix = createLaplacianOperator(number_points);
}

Eigen::MatrixXd Discretization::createLaplacianOperator(const int number_points)
{
    Eigen::VectorXd diagonal_Laplace = -2 * Eigen::VectorXd::Ones(number_points);

    Eigen::VectorXd upper_diagonal_Laplace = Eigen::VectorXd::Ones(number_points - 1);

    Eigen::VectorXd lower_diagonal_Laplace = Eigen::VectorXd::Ones(number_points - 1);

    Eigen::MatrixXd laplacian_matrix = Eigen::MatrixXd::Zero(number_points, number_points);

    laplacian_matrix.diagonal() = diagonal_Laplace;
    laplacian_matrix.diagonal(1) = upper_diagonal_Laplace;
    laplacian_matrix.diagonal(-1) = lower_diagonal_Laplace;

    return laplacian_matrix;
}