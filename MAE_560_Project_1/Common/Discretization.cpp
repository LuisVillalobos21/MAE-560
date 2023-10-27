#include "Discretization.hpp"

Eigen::MatrixXd createLaplacianMatrix(const int number_points)
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

Eigen::MatrixXd createDifferentiationMatrix(const int number_points)
{
    Eigen::VectorXd diagonal = Eigen::VectorXd::Zero(number_points);

    Eigen::VectorXd upper_diagonal = Eigen::VectorXd::Ones(number_points - 1);

    Eigen::VectorXd lower_diagonal = - Eigen::VectorXd::Ones(number_points - 1);

    Eigen::MatrixXd differentiaion_matrix = Eigen::MatrixXd::Zero(number_points, number_points);

    differentiaion_matrix.diagonal() = diagonal;
    differentiaion_matrix.diagonal(1) = upper_diagonal;
    differentiaion_matrix.diagonal(-1) = lower_diagonal;

    return differentiaion_matrix;
}

///////// EXPLICIT EULER HEAT EQN
Discretization::Discretization(const int number_points)
{
    laplacian_matrix = createLaplacianMatrix(number_points);
}

///////// CRANK NICOLSON HEAT EQN
CNDiscretization::CNDiscretization(const int number_points, const Parameters& params, const Mesh& mesh)
{
    laplacian_matrix = createLaplacianMatrix(number_points);
    cn_matrix = createCNMatrix(params, mesh);
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

///////// LAX WENDROFF WAVE EQN
LWDiscretization::LWDiscretization(const int number_points)
{
    laplacian_matrix = createLaplacianMatrix(number_points);
    differentiaion_matrix = createDifferentiationMatrix(number_points);
}