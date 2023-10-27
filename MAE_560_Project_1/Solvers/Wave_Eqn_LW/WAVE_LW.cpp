#include "ProjectIncludes.hpp"

// LAX WENDROFF FOR WAVE EQN

void writeToFile(std::ofstream& file, const Eigen::VectorXd& u)
{
    file << u.transpose() << std::endl;
}

Eigen::VectorXd computeInitialCondition(const Eigen::VectorXd& x)
{
    Eigen::VectorXd u = Eigen::VectorXd::Zero(x.size());

    for (int i = 0; i < x.size(); ++i)
    {
        u(i) = std::exp(-(x(i) - 0.5) * (x(i) - 0.5) / 0.01);
    }

    return u;
}

int main()
{
    double domainStart = 0.0;
    double domainEnd = 1.0;
    double dx = 0.01;

    PeriodicMesh mesh(domainStart, domainEnd, dx);

    WaveParameters params;
    double c = 1;
    double CFL = .4;
    double dt{ (CFL * dx) / c };

    params.dt = dt;
    params.c = c;

    SolutionVector soln(mesh.number_points);
    Eigen::VectorXd initialCondition = computeInitialCondition(mesh.x);
    soln.u = initialCondition;

    PeriodicBoundaryConditions boundaries(params, mesh, soln);

    LWDiscretization discrete(mesh.number_points);

    LWTimeIntegration LWupdate;

    std::ofstream file("AB2_initial_data_CFL_4.csv");

    //writeToFile(file, soln.u);

    int number_time_steps = 2;
    int step = 0;
    auto start_time = std::chrono::high_resolution_clock::now();
    for(int i = 1; i < number_time_steps; ++i)
    {
        boundaries.calcPeriodicValue(soln,mesh);
        LWupdate.laxWendroff(params, mesh, soln, boundaries, discrete);

        //std::cout << "Time step: " << i << '\n';
        writeToFile(file, soln.u);

        ++step;
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_time - start_time;

    file.close();

    std::cout << "Solution u is: " << soln.u << '\n';
    std::cout << "Mesh is: " << mesh.x << '\n';
    std::cout << "Simulation ended at time step: " << step << '\n';
    std::cout << "Elapsed time: " << elapsed_seconds.count() << "s\n";

    return 0;
}
