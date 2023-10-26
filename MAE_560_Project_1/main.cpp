#include "ProjectIncludes.hpp"

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

    //std::ofstream file("u_values.csv");

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
