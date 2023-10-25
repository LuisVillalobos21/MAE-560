#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <fstream>

struct SimulationParameters
{
    double dx = 0.01;
    double domain_start = 0.0;
    double domain_end = 1.0;
    int number_points = static_cast<int>((domain_end + dx - domain_start) / dx);
    double alpha = 5e-3;
    double dt = .009;
    double a = 0.4;
    double b = -0.1;
    double c = 0.1;
    double u_bc_L = 1;

    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(number_points, domain_start, domain_end);
    Eigen::VectorXd u = Eigen::VectorXd::Zero(number_points);
    Eigen::VectorXd u_new = Eigen::VectorXd::Zero(number_points);
    Eigen::VectorXd u_bc = Eigen::VectorXd::Zero(number_points);

    Eigen::MatrixXd createBCMatrix()
    {
        Eigen::MatrixXd bc_matrix = Eigen::MatrixXd::Zero(2, 2);
        bc_matrix(0, 0) = 2;
        bc_matrix(0, 1) = -1;
        bc_matrix(1, 0) = -(2 * a * dx) / b;
        bc_matrix(1, 1) = 1;
        return bc_matrix;
    }
};

// Laplacian operator creation
Eigen::MatrixXd createLaplacianOperator(const int number_points)
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

// Calculate the robin u_bc value
double calcRobinValue(const SimulationParameters& params, const Eigen::MatrixXd bc_matrix)
{
    Eigen::VectorXd bc_rhs = Eigen::VectorXd::Zero(2);
    bc_rhs(0) = params.u(0);
    bc_rhs(1) = params.u(0) - (2 * params.c * params.dx) / params.b;

    Eigen::VectorXd f_bc = bc_matrix.partialPivLu().solve(bc_rhs);

    return f_bc(0);
}

// Explict euler integration
Eigen::MatrixXd explicitEuler(SimulationParameters& params, const Eigen::MatrixXd laplacian_matrix, const Eigen::MatrixXd bc_matrix)
{
    params.u_bc(0) = calcRobinValue(params, bc_matrix);
    params.u_bc(params.u_bc.size() - 1) = params.u_bc_L;

    Eigen::MatrixXd u_new = params.u + ((params.alpha * params.dt) / (params.dx * params.dx)) * (laplacian_matrix * params.u + params.u_bc);

    return u_new;
}

// Data out function
void writeToFile(std::ofstream& file, const Eigen::VectorXd& u)
{
    file << u.transpose() << std::endl;
}


int main()
{
    SimulationParameters params;

    Eigen::MatrixXd bc_matrix = params.createBCMatrix();

    std::cout << "Inital u is: " << params.u << std::endl;

    std::ofstream file("u_values.csv");

    writeToFile(file, params.u);

    Eigen::MatrixXd laplacian_matrix = createLaplacianOperator(params.number_points);

    // Time loop
    int number_time_steps{ 15000 };
    for (int i = 1; i <= number_time_steps; ++i) 
    {
        params.u = explicitEuler(params,laplacian_matrix,bc_matrix);

        writeToFile(file, params.u);

        std::cout << "Completed time step " << i << '\n';
    }

    
    std::cout << "Solution Vector is:\n" << params.u << std::endl;

    file.close(); // Close the file

    return 0;
}
