#include "ProjectIncludes.hpp"

// EXPLICIT EULER SCHEME FOR HEAT EQN

void writeToFile(std::ofstream& file, const Eigen::VectorXd& u)
{
    file << u.transpose() << std::endl;
}

Eigen::VectorXd analyticalSolution(const Eigen::VectorXd& x_values) {
    Eigen::VectorXd result = Eigen::VectorXd::Constant(x_values.size(), 0.4);
    result += 0.6 * x_values;
    return result;
}

double computeRelativeError(const Eigen::VectorXd& computed, const Eigen::VectorXd& analytical)
{
    Eigen::VectorXd diff = computed - analytical;
    double norm_diff = diff.norm();
    double norm_analytical = analytical.norm();

    return norm_diff / norm_analytical;
}

int main()
{
    double domainStart = 0.0;
    double domainEnd = 1.0;
    double dx = 0.01;

    Mesh mesh(domainStart, domainEnd, dx);

    Parameters params;
    double alpha{.005};
    double dt{.01};
    double a = 0.4;
    double b = -0.1;
    double c = 0.1;

    params.alpha = alpha;
    params.dt = dt;
    params.a = a;
    params.b = b;
    params.c = c;

    SolutionVector soln(mesh.number_points);

    BoundaryConditions boundaries(mesh.number_points, params, mesh, soln);

    Discretization discrete(mesh.number_points);

    TimeIntegration euler_update;

    std::ofstream file("heat_explicit.csv");

    //writeToFile(file, soln.u);

    Eigen::VectorXd u_anal = analyticalSolution(mesh.x);
    double error = computeRelativeError(soln.u, u_anal);

    int i = 1;
    auto start_time = std::chrono::high_resolution_clock::now();
    while (error > 1e-3)
    {
        boundaries.calcRobinValueExp(params,mesh, soln);
        euler_update.explicitEuler(params, mesh, soln, boundaries, discrete);

        error = computeRelativeError(soln.u, u_anal);
        //std::cout << "Time step: " << i << ", Error: " << error << '\n';
        //writeToFile(file, soln.u);

        ++i;
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_time - start_time;

    Eigen::VectorXd uplot = Eigen::VectorXd::Zero(mesh.number_points + 2);
    uplot(0) = ((params.c - (params.b / mesh.dx) * soln.u(0)) / (params.a - (params.b / mesh.dx)));
    uplot.segment(1, mesh.number_points) = soln.u;
    uplot(uplot.size() - 1) = 1;

    writeToFile(file, uplot);
    file.close();

    std::cout << "Mesh is " << mesh.x << '\n';
    //std::cout << "Solution u is: " << soln.u << '\n';
    std::cout << "Solution u is: " << uplot << '\n';
    std::cout << "Simulation ended at time step: " << i << ", with error: " << error << '\n';
    std::cout << "Elapsed time: " << elapsed_seconds.count() << "s\n";
    
    return 0;
}
