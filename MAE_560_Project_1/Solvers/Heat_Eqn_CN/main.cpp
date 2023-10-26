#include "ProjectIncludes.hpp"

 // CRANK NICOLSON SCHEME FOR HEAT EQN

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
    double alpha{ .005 };
    double dt{ .1 };
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

    CNDiscretization discrete(mesh.number_points,params,mesh);

    FormLinearSolveEqn rhs(mesh.number_points, params, mesh, soln, boundaries, discrete);

    CNTimeIntegration cn_update;

    //std::ofstream file("u_values.csv");

    //writeToFile(file, soln.u);

    Eigen::VectorXd u_anal = analyticalSolution(mesh.x);
    double error = computeRelativeError(soln.u, u_anal);
    int i = 1;

    while (error > 1e-2)
    {
        boundaries.calcRobinValue(params, mesh, soln);
        rhs.calcLinearSolveRHS(params, mesh, soln, boundaries, discrete);
        cn_update.thomasAlgorithm(discrete, rhs, soln);

        error = computeRelativeError(soln.u, u_anal);

        std::cout << "Time step: " << i << ", Error: " << error << '\n';
        //writeToFile(file, soln.u);

        ++i;
    }

    //file.close();

    std::cout << "Solution u is: " << soln.u << '\n';
    std::cout << "Simulation ended at time step: " << i << ", with error: " << error << '\n';

    return 0;
}