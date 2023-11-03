#include "ProjectIncludes.hpp"

// AB2 FOR WAVE EQN

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
    double CFL = 0.4;
    double dt{ (CFL * dx) / c };

    params.dt = dt;
    params.c = c;

    SolutionVector soln(mesh.number_points);
    Eigen::VectorXd initialCondition = computeInitialCondition(mesh.x);
    soln.u = readCSVtoEigenVector("AB2_initial_data_CFL_4.csv");
    soln.u_old = initialCondition;

    AB2PeriodicBoundaryConditions boundaries(params, mesh, soln);

    LWDiscretization discrete(mesh.number_points);

    AB2TimeIntegration AB2update;

    std::ofstream file("AB2_CFL_4.csv");

    //writeToFile(file, soln.u);

    int number_time_steps = 1001;
    int step = 0;
    auto start_time = std::chrono::high_resolution_clock::now();
    for (int i = 1; i < number_time_steps; ++i)
    {
        //std::cout << "Current u is " << soln.u << '\n';
        //std::cout << "Old u is " << soln.u_old << '\n';
        boundaries.calcPeriodicValueAB2(soln, mesh);
        AB2update.AB2(params, mesh, soln, boundaries, discrete);
        //std::cout << "Updated u is " << soln.u << '\n';
        //std::cout << "New old u is " << soln.u_old << '\n';

        //std::cout << "Time step: " << i << '\n';
        writeToFile(file, soln.u);

        ++step;
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_time - start_time;
    //writeToFile(file, soln.u);

    file.close();

    std::cout << "Solution u is: " << soln.u << '\n';
    std::cout << "Simulation ended at time step: " << step << '\n';
    std::cout << "Elapsed time: " << elapsed_seconds.count() << "s\n";

    return 0;
}
