#include "Discretization.hpp"

///////// EXPLICIT EULER HEAT EQN

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

///////// CRANK NICOLSON HEAT EQN

CNDiscretization::CNDiscretization(const int number_points, const Parameters& params, const Mesh& mesh) 
{
    laplacian_matrix = createLaplacianOperator(number_points);
    cn_matrix = createCNMatrix(params, mesh);
}

Eigen::MatrixXd CNDiscretization::createLaplacianOperator(const int number_points)
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

Eigen::MatrixXd CNDiscretization::createCNMatrix(const Parameters& params, const Mesh& mesh)
{
    Eigen::VectorXd diagonal = (1 + (params.alpha * params.dt) / (mesh.dx * mesh.dx)) *
                               Eigen::VectorXd::Ones(mesh.number_points);
    
    Eigen::VectorXd upper_diagonal = ((-params.alpha * params.dt) / (2 * mesh.dx * mesh.dx)) *
                                     Eigen::VectorXd::Ones(mesh.number_points - 1);

    Eigen::VectorXd lower_diagonal = upper_diagonal;

    Eigen::MatrixXd cn_lhs_matrix = Eigen::MatrixXd::Zero(mesh.number_points, mesh.number_points);

    cn_lhs_matrix.diagonal() = diagonal;
    cn_lhs_matrix.diagonal(1) = upper_diagonal;
    cn_lhs_matrix.diagonal(-1) = lower_diagonal;

    return cn_lhs_matrix;
}
