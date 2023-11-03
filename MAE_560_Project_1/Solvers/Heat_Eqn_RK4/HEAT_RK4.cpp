#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <cmath>
const double PI = 3.14159265358979323846;

// HEAT EQUATION WITH RK4 AND SOURCE TERM

struct MeshStruct
{
    double dx;
    double domain_start;
    double domain_end;
    int number_points;
    Eigen::VectorXd x;

    MeshStruct(double start, double end, double step)
        : domain_start(start), domain_end(end), dx(step)
    {
        number_points = static_cast<int>((domain_end + dx - domain_start) / dx) - 2;
        x = Eigen::VectorXd::LinSpaced(number_points, domain_start, domain_end - dx);
    }

};

struct ParameterStruct
{
    double alpha;
    double dt;
    double k;
    double omega;
    double a;
    double b;
    double c;
    double d;
    double P;
    double Q;

    ParameterStruct(double alpha, double dt, double k, double omega, double a, double b, double c, double d, const MeshStruct& mesh)
        : alpha(alpha), dt(dt), k(k), omega(omega), a(a), b(b), c(c), d(d),
        P((-alpha * dt) / (2 * mesh.dx * mesh.dx)), Q(1 + (alpha * dt) / (mesh.dx * mesh.dx)) {}
};

struct SolutionStruct
{
    Eigen::VectorXd u;

    SolutionStruct(const MeshStruct& mesh) : u(Eigen::VectorXd::Zero(mesh.number_points)) 
    {

    }
};

struct OperatorStruct
{
    Eigen::MatrixXd laplacian_matrix;
    Eigen::MatrixXd differentiation_matrix;

    OperatorStruct(const MeshStruct mesh, const ParameterStruct params)
    {
        laplacian_matrix = createLaplacianMatrix(mesh, params);
        differentiation_matrix = createDifferentiationMatrix(mesh, params);
    }

    Eigen::MatrixXd createLaplacianMatrix(const MeshStruct mesh, const ParameterStruct params)
    {
        Eigen::VectorXd diagonal_Laplace = (1 / (mesh.dx * mesh.dx)) * -2 * Eigen::VectorXd::Ones(mesh.number_points);

        Eigen::VectorXd upper_diagonal_Laplace = (1 / (mesh.dx * mesh.dx)) * Eigen::VectorXd::Ones(mesh.number_points - 1);

        Eigen::VectorXd lower_diagonal_Laplace = (1 / (mesh.dx * mesh.dx)) * Eigen::VectorXd::Ones(mesh.number_points - 1);

        Eigen::MatrixXd laplacian_matrix = Eigen::MatrixXd::Zero(mesh.number_points, mesh.number_points);

        laplacian_matrix.diagonal() = diagonal_Laplace;
        laplacian_matrix.diagonal(1) = upper_diagonal_Laplace;
        laplacian_matrix.diagonal(-1) = lower_diagonal_Laplace;

        return laplacian_matrix;
    }

    Eigen::MatrixXd createDifferentiationMatrix(const MeshStruct mesh, const ParameterStruct params)
    {
        Eigen::VectorXd upper_diagonal = (0.5 / mesh.dx) * Eigen::VectorXd::Ones(mesh.number_points - 1);

        Eigen::VectorXd lower_diagonal = (0.5 / mesh.dx) * -Eigen::VectorXd::Ones(mesh.number_points - 1);

        Eigen::MatrixXd differentiaion_matrix = Eigen::MatrixXd::Zero(mesh.number_points, mesh.number_points);

        differentiaion_matrix.diagonal(1) = upper_diagonal;
        differentiaion_matrix.diagonal(-1) = lower_diagonal;

        return differentiaion_matrix;
    }
};

struct BCStruct
{
    Eigen::VectorXd u_bc;

    BCStruct(const ParameterStruct& params, const MeshStruct mesh, const SolutionStruct soln)
        : u_bc(Eigen::VectorXd::Zero(mesh.number_points))
    {
        updateBC(params, mesh, soln);
    }

    void updateBC(const ParameterStruct& params, const MeshStruct mesh, const SolutionStruct& soln)
    {
        u_bc(0) = (params.alpha / (mesh.dx * mesh.dx)) * ((params.c - params.b * (soln.u(0) / mesh.dx)) / (params.a - (params.b / mesh.dx)));

        u_bc(u_bc.size() - 1) = (params.alpha / (mesh.dx * mesh.dx)) * (soln.u(soln.u.size() - 1) + mesh.dx * params.d);
    }
};

struct ksourceStruct
{
    Eigen::VectorXd ksource;
    const MeshStruct& mesh;
    const ParameterStruct& params;

    ksourceStruct(const MeshStruct& mesh, const ParameterStruct& params)
        : ksource(Eigen::VectorXd::Zero(mesh.number_points)), mesh(mesh), params(params){}

    void calckSource(double t)
    {
        ksource = params.k * std::sin(params.omega * t) * (-(mesh.x.array() - 0.5) * (mesh.x.array() - 0.5) / 0.01).exp();
    }
};

struct RK4Struct
{
    Eigen::VectorXd g;
    Eigen::VectorXd k1;
    Eigen::VectorXd k2;
    Eigen::VectorXd k3;
    Eigen::VectorXd k4;

    const MeshStruct& mesh;
    const ParameterStruct& params;
    SolutionStruct& soln;
    ksourceStruct& ksrc;
    const OperatorStruct& op;
    BCStruct& BCs;

    RK4Struct(const MeshStruct& mesh, const ParameterStruct& params, SolutionStruct& soln, ksourceStruct& ksrc,
        const OperatorStruct& op, BCStruct& BCs)
        : mesh(mesh), params(params), soln(soln), ksrc(ksrc), op(op), BCs(BCs),
        k1(Eigen::VectorXd::Zero(mesh.number_points)),
        k2(Eigen::VectorXd::Zero(mesh.number_points)),
        k3(Eigen::VectorXd::Zero(mesh.number_points)),
        k4(Eigen::VectorXd::Zero(mesh.number_points)),
        g(Eigen::VectorXd::Zero(mesh.number_points))
    {

    }

    Eigen::VectorXd gUpdate(double time, Eigen::VectorXd soln)
    {
        ksrc.calckSource(time);
        g = params.alpha * (op.laplacian_matrix * soln) + ksrc.ksource + BCs.u_bc;
        return g;
    }

    void calck1(double time, Eigen::VectorXd soln)
    {
        k1 = gUpdate(time, soln);
        //std::cout << "Source vector is" << ksrc.ksource << '\n';
        //std::cout << "BC vector is" << BCs.u_bc << '\n';
        //std::cout << "Gupdate is " << params.alpha * (op.laplacian_matrix * soln) + ksrc.ksource + BCs.u_bc << '\n';
    }

    void calck2(double time, Eigen::VectorXd soln)
    {
        k2 = gUpdate(time + 0.5 * params.dt, soln + 0.5 * params.dt * k1);
    }

    void calck3(double time, Eigen::VectorXd soln)
    {
        k3 = gUpdate(time + 0.5 * params.dt, soln + 0.5 * params.dt * k2);
    }

    void calck4(double time, Eigen::VectorXd soln)
    {
        k4 = gUpdate(time + params.dt, soln + params.dt * k3);
    }

    void timeMarchRK4(double time)
    {
        BCs.updateBC(params, mesh, soln);

        calck1(time, soln.u);
        //std::cout << "K1 is " << k1 << '\n';
        calck2(time, soln.u);
        //std::cout << "K2 is " << k2 << '\n';
        calck3(time, soln.u);
        //std::cout << "K3 is " << k3 << '\n';
        calck4(time, soln.u);
        //std::cout << "K4 is " << k4 << '\n';

        //std::cout << "Update to U is " << (1.0 / 6.0) * params.dt * (k1 + k2 + k3 + k4) << '\n';
        //std::cout << "Old u is " << soln.u << '\n';
        soln.u = soln.u + (1.0 / 6.0) * params.dt * (k1 + k2 + k3 + k4);
        //std::cout << "New u is  " << soln.u << '\n';
    }
};

void writeToFile(std::ofstream& file, const Eigen::VectorXd& u)
{
    file << u.transpose() << std::endl;
}

int main()
{
    MeshStruct mesh(0, 1, .01);
    ParameterStruct params(.005, .01, 0.5, 0.2 * PI, 0.4, -0.1, 0.1, 0.6, mesh);
    SolutionStruct soln(mesh);
    OperatorStruct op(mesh, params);
    BCStruct BCs(params, mesh, soln);
    ksourceStruct src(mesh, params);
    RK4Struct rk(mesh, params, soln, src, op, BCs);

    std::cout << "Intial u is " << '\n' << soln.u << '\n';
    std::ofstream file("heatrk4_t196.csv");
    //writeToFile(file, soln.u);

    double time = 0;
    for (int i = 1; i < 19601; ++i)
    {
        time += params.dt;
        rk.timeMarchRK4(time);
        //writeToFile(file, soln.u);
    }

    Eigen::VectorXd uplot = Eigen::VectorXd::Zero(mesh.number_points + 2);
    uplot(0) = ((params.c - (params.b / mesh.dx) * soln.u(0)) / (params.a - (params.b / mesh.dx)));
    uplot.segment(1, mesh.number_points) = soln.u;
    uplot(uplot.size() - 1) = params.d * mesh.dx + soln.u(soln.u.size() - 1);

    //std::cout << "Solution u is " << '\n' << soln.u << '\n';

    writeToFile(file, uplot);
    file.close();

	return 0;
}