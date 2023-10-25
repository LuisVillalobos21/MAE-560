#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <fstream>

struct Mesh
{
    double dx;
    double domain_start;
    double domain_end;
    int number_points;
    Eigen::VectorXd x;

    Mesh(double start, double end, double step) 
    {
        domain_start = start;
        domain_end = end;
        dx = step;
        number_points = static_cast<int>((domain_end + dx - domain_start) / dx);
        x = Eigen::VectorXd::LinSpaced(number_points, domain_start, domain_end);
    }

};

struct Parameters
{
    double alpha;
    double dt;
    double a;
    double b;
    double c;
};

struct SolutionVector
{
    Eigen::VectorXd u;
    Eigen::VectorXd u_new;
    SolutionVector(int number_points)
    {
        u = Eigen::VectorXd::Zero(number_points);
        u_new = Eigen::VectorXd::Zero(number_points);
    }
};

struct BoundaryConditions
{
    Eigen::VectorXd u_bc;
    Eigen::MatrixXd robin_matrix;
    Eigen::VectorXd bc_rhs;
    double u_bc_L = 1;

    BoundaryConditions(int number_points, const Parameters& params, const Mesh& mesh,
        const SolutionVector& soln)
    {
        u_bc = Eigen::VectorXd::Zero(number_points);
        bc_rhs = Eigen::VectorXd::Zero(2);
        robin_matrix = createRobinLeftBCMatrix(params,mesh);
        calcRobinValue(params, mesh, soln);
        u_bc(u_bc.size() - 1) = u_bc_L;
    }

    Eigen::MatrixXd createRobinLeftBCMatrix(const Parameters& params, const Mesh& mesh)
    {
        Eigen::MatrixXd bc_matrix = Eigen::MatrixXd::Zero(2, 2);
        bc_matrix(0, 0) = 2;
        bc_matrix(0, 1) = -1;
        bc_matrix(1, 0) = -(2 * params.a * mesh.dx) / params.b;
        bc_matrix(1, 1) = 1;
        return bc_matrix;
    }

    void calcRobinValue(const Parameters& params, const Mesh& mesh, const SolutionVector& soln)
    {
        bc_rhs(0) = soln.u(0);
        bc_rhs(1) = soln.u(0) - (2 * params.c * mesh.dx) / params.b;

        Eigen::VectorXd f_bc = robin_matrix.partialPivLu().solve(bc_rhs);
        u_bc(0) = f_bc(0);
    }
};

struct Discretization
{
    Eigen::MatrixXd laplacian_matrix;

    Discretization(const int number_points)
    {
        laplacian_matrix = createLaplacianOperator(number_points);
    }
        
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
};

struct TimeIntegration
{
    void explicitEuler(const Parameters& params, const Mesh& mesh, SolutionVector& soln, 
        const BoundaryConditions& boundaries, const Discretization& discrete)
    {
        soln.u_new = soln.u + ((params.alpha * params.dt) / (mesh.dx * mesh.dx)) * (discrete.laplacian_matrix * soln.u + boundaries.u_bc);
        soln.u = soln.u_new;
    }
};

void writeToFile(std::ofstream& file, const Eigen::VectorXd& u)
{
    file << u.transpose() << std::endl;
}

int main()
{
    double domainStart = 0.0;
    double domainEnd = 1.0;
    double dx = 0.01;

    Mesh mesh(domainStart, domainEnd, dx);

    Parameters params;
    double alpha{.005};
    double dt{.0099};
    double a = 0.4;
    double b = -0.1;
    double c = 0.1;

    params.alpha = alpha;
    params.dt = dt;
    params.a = a;
    params.b = b;
    params.c = c;

    SolutionVector soln(mesh.number_points);

    std::ofstream file("u_values.csv");

    //writeToFile(file, soln.u);

    BoundaryConditions boundaries(mesh.number_points, params, mesh, soln);

    Discretization discrete(mesh.number_points);

    TimeIntegration euler_update;

    int number_time_steps = 15000;

    for (int i = 1; i < number_time_steps; ++i)
    {
        boundaries.calcRobinValue(params,mesh, soln);
        euler_update.explicitEuler(params, mesh, soln, boundaries, discrete);
        std::cout << "Time step: " << i <<'\n';
        //writeToFile(file, soln.u);
    }

    //file.close();

    std::cout << soln.u << '\n';
   
    return 0;
}
