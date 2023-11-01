#include <iostream>
#include <fstream>
#include <Eigen/Dense>

// ADVECTION DIFFUSION EQUATION

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
        number_points = static_cast<int>((domain_end + dx - domain_start) / dx) - 1;
        x = Eigen::VectorXd::LinSpaced(number_points, domain_start, domain_end - dx);
    }

};

struct ParameterStruct
{
    double alpha;
    double c;
    double dt;
    double P;
    double Q;

    ParameterStruct(double alpha, double c, double dt,const MeshStruct& mesh)
        : alpha(alpha), c(c), dt(dt), P((-alpha * dt) / (2 * mesh.dx * mesh.dx)), Q(1 + (alpha * dt) / (mesh.dx * mesh.dx)){}
};

struct SolutionStruct
{
    Eigen::VectorXd u;
    Eigen::VectorXd u_new;
    Eigen::VectorXd u_old;
    SolutionStruct(const MeshStruct& mesh)
    {
        Eigen::VectorXd u = Eigen::VectorXd::Zero(mesh.number_points);
        Eigen::VectorXd u_new = Eigen::VectorXd::Zero(mesh.number_points);
        Eigen::VectorXd u_old = Eigen::VectorXd::Zero(mesh.number_points);
    }
};

struct BCStruct
{
    Eigen::VectorXd u_bc_diff_AB2;
    Eigen::VectorXd u_bc_diff_exp;
    Eigen::VectorXd u_bc_laplacian;

    BCStruct(const ParameterStruct& params, const MeshStruct mesh, const SolutionStruct soln)
        : u_bc_diff_AB2(Eigen::VectorXd::Zero(mesh.x.size())), u_bc_laplacian(Eigen::VectorXd::Zero(mesh.x.size())),
          u_bc_diff_exp(Eigen::VectorXd::Zero(mesh.x.size()))
    {
        updateBC_diff_exp(params, mesh, soln);
        updateBC_lap(params, mesh, soln);
    }

    void updateBC_diff_AB2(const ParameterStruct& params, const MeshStruct mesh, const SolutionStruct& soln)
    {
        u_bc_diff_AB2(0) = (1.5 * params.c * params.dt) * ((0.5 * soln.u(soln.u.size() - 1)) / mesh.dx) - 
            (0.5 * params.c * params.dt) * (0.5 * soln.u_old(soln.u_old.size() - 1) / mesh.dx);

        u_bc_diff_AB2(soln.u.size() - 1) = (-1.5 * params.c * params.dt) * ((0.5 * soln.u(0)) / mesh.dx) + 
            (0.5 * params.c * params.dt) * (0.5 * soln.u_old(0) / mesh.dx);
    }
    void updateBC_diff_exp(const ParameterStruct& params, const MeshStruct mesh, const SolutionStruct& soln)
    {
        u_bc_diff_exp(0) = (params.c * params.dt) * ((0.5 * soln.u(soln.u.size() - 1)) / mesh.dx);
        u_bc_diff_exp(soln.u.size() - 1) = (-params.c * params.dt) * ((0.5 * soln.u(0)) / mesh.dx);
    }
    void updateBC_lap(const ParameterStruct& params, const MeshStruct mesh, const SolutionStruct& soln)
    {
        u_bc_laplacian(0) = (0.5 * params.dt) * ((params.alpha * soln.u(soln.u.size() - 1)) / (mesh.dx * mesh.dx));
        u_bc_laplacian(soln.u.size() - 1) = (0.5 * params.dt) * ((params.alpha * soln.u(0)) / (mesh.dx * mesh.dx));
    }
};

struct CN_MatrixStruct
{
    Eigen::MatrixXd CN_Matrix;
    CN_MatrixStruct(const ParameterStruct& params, const MeshStruct& mesh)
        :CN_Matrix(Eigen::MatrixXd::Zero(mesh.number_points - 1, mesh.number_points - 1))
    {
        constructMatrix(params, mesh);
    }

    void constructMatrix(const ParameterStruct& params, const MeshStruct& mesh)
    {
        Eigen::VectorXd diagonal = (1 + (params.alpha * params.dt) / (mesh.dx * mesh.dx)) * Eigen::VectorXd::Ones(mesh.number_points -1);

        Eigen::VectorXd upper_diagonal = ((-params.alpha * params.dt) / (2 * mesh.dx * mesh.dx)) * Eigen::VectorXd::Ones(mesh.number_points - 2);

        Eigen::VectorXd lower_diagonal = upper_diagonal;

        CN_Matrix.diagonal() = diagonal;
        CN_Matrix.diagonal(1) = upper_diagonal;
        CN_Matrix.diagonal(-1) = lower_diagonal;
    }

};

struct OperatorStruct
{
    Eigen::MatrixXd laplacian_matrix;
    Eigen::MatrixXd differentiation_matrix;

    OperatorStruct(const MeshStruct mesh, const ParameterStruct params)
    {
        laplacian_matrix = createLaplacianMatrix(mesh,params);
        differentiation_matrix = createDifferentiationMatrix(mesh,params);
    }

    Eigen::MatrixXd createLaplacianMatrix(const MeshStruct mesh, const ParameterStruct params)
    {
        Eigen::VectorXd diagonal_Laplace = (params.alpha / (mesh.dx * mesh.dx)) * - 2 * Eigen::VectorXd::Ones(mesh.number_points);

        Eigen::VectorXd upper_diagonal_Laplace = (params.alpha / (mesh.dx * mesh.dx)) * Eigen::VectorXd::Ones(mesh.number_points - 1);

        Eigen::VectorXd lower_diagonal_Laplace = (params.alpha / (mesh.dx * mesh.dx)) * Eigen::VectorXd::Ones(mesh.number_points - 1);

        Eigen::MatrixXd laplacian_matrix = Eigen::MatrixXd::Zero(mesh.number_points, mesh.number_points);

        laplacian_matrix.diagonal() = diagonal_Laplace;
        laplacian_matrix.diagonal(1) = upper_diagonal_Laplace;
        laplacian_matrix.diagonal(-1) = lower_diagonal_Laplace;

        return laplacian_matrix;
    }

    Eigen::MatrixXd createDifferentiationMatrix(const MeshStruct mesh, const ParameterStruct params)
    {
        Eigen::VectorXd upper_diagonal = (0.5 / mesh.dx) * Eigen::VectorXd::Ones(mesh.number_points - 1);

        Eigen::VectorXd lower_diagonal = (0.5 / mesh.dx)  * -Eigen::VectorXd::Ones(mesh.number_points - 1);

        Eigen::MatrixXd differentiaion_matrix = Eigen::MatrixXd::Zero(mesh.number_points, mesh.number_points);

        differentiaion_matrix.diagonal(1) = upper_diagonal;
        differentiaion_matrix.diagonal(-1) = lower_diagonal;

        return differentiaion_matrix;
    }
};

struct RHSstruct
{
    Eigen::VectorXd d_AB2;
    Eigen::VectorXd d_exp;
    Eigen::VectorXd q;

    RHSstruct(const ParameterStruct& params, const MeshStruct mesh, const SolutionStruct& soln, const OperatorStruct op, const BCStruct boundaries)
        : d_AB2(Eigen::VectorXd::Zero(mesh.number_points)), d_exp(Eigen::VectorXd::Zero(mesh.number_points)), q(Eigen::VectorXd::Zero(mesh.number_points - 1))
    {
        d_exp_update(params, mesh, soln, op, boundaries);
        q_update(params);
    }

    void d_AB2_update(const ParameterStruct& params, const MeshStruct mesh, const SolutionStruct& soln, const OperatorStruct op, const BCStruct boundaries)
    {
        d_AB2 = soln.u + (-1.5 * params.c * params.dt) * (op.differentiation_matrix * soln.u) + (0.5 * params.c * params.dt) * (op.differentiation_matrix * soln.u_old) +
            (0.5 * params.dt) * (op.laplacian_matrix * soln.u) + boundaries.u_bc_diff_AB2 + boundaries.u_bc_laplacian;
    }

    void d_exp_update(const ParameterStruct& params, const MeshStruct mesh, const SolutionStruct& soln, const OperatorStruct op, const BCStruct boundaries)
    {
        d_exp = soln.u + (-params.c * params.dt) * (op.differentiation_matrix * soln.u) +
            (0.5 * params.dt) * (op.laplacian_matrix * soln.u) + boundaries.u_bc_diff_exp + boundaries.u_bc_laplacian;
    }

    void q_update(const ParameterStruct& params)
    {
        q(0) = -params.P;
        q(q.size() - 1) = -params.P;
    }

};

Eigen::VectorXd thomasAlgorithm(const Eigen::MatrixXd& matrix, const Eigen::VectorXd& rhs)
{
    int n = rhs.size();

    Eigen::VectorXd a(n), b_prime(n), c(n), d = rhs, soln(n);

    for (int i = 0; i < n; i++) {
        b_prime(i) = matrix(i, i);
        if (i > 0) a(i) = matrix(i, i - 1);
        if (i < n - 1) c(i) = matrix(i, i + 1);
    }

    for (int i = 1; i < n; i++) {
        double m = a(i) / b_prime(i - 1);
        b_prime(i) = b_prime(i) - m * c(i - 1);
        d(i) = d(i) - m * d(i - 1);
    }

    soln(n - 1) = d(n - 1) / b_prime(n - 1);
    for (int i = n - 2; i >= 0; i--) {
        soln(i) = (d(i) - c(i) * soln(i + 1)) / b_prime(i);
    }

    return soln;
}

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
    MeshStruct mesh(0, 1, .01);

    ParameterStruct params(.005, 1, .0075, mesh);

    SolutionStruct soln(mesh);

    Eigen::VectorXd initialCondition = computeInitialCondition(mesh.x);
    soln.u = initialCondition;
    std::ofstream file("adv_diff.csv");
    writeToFile(file, soln.u);
    //std::cout << "U intial is " << '\n' << soln.u << '\n';

    OperatorStruct op(mesh, params);

    BCStruct boundaries(params,mesh,soln);

    CN_MatrixStruct matrix(params, mesh);

    RHSstruct rhs(params, mesh, soln, op, boundaries);
    
    Eigen::VectorXd X2;
    for (int i = 1; i < 2; ++i)
    {
        boundaries.updateBC_diff_exp(params, mesh, soln);
        boundaries.updateBC_lap(params, mesh, soln);
        rhs.d_exp_update(params, mesh, soln, op, boundaries);

        soln.u_old = soln.u;
        soln.u.head(soln.u.size() - 1) = thomasAlgorithm(matrix.CN_Matrix, rhs.d_exp.head(rhs.d_exp.size() - 1));
        X2 = thomasAlgorithm(matrix.CN_Matrix, rhs.q);

        soln.u(soln.u.size() - 1) = (rhs.d_exp(rhs.d_exp.size() - 1) - params.P * soln.u(0) - params.P * soln.u(soln.u.size() - 2)) /
            (params.Q + params.P * X2(0) + params.P * X2(X2.size() - 1));

        soln.u.head(soln.u.size() - 1) = soln.u.head(soln.u.size() - 1) + X2 * soln.u(soln.u.size() - 1);
        
        writeToFile(file, soln.u);
    }
    
    //std::cout << "U old from explicit " << '\n' << soln.u_old << '\n';
    //std::cout << "U from explicit " << '\n' << soln.u << '\n';

    for (int i = 1; i < 20000; ++i)
    {
        boundaries.updateBC_diff_AB2(params, mesh, soln);
        boundaries.updateBC_lap(params, mesh, soln);
        //if (i == 999)
        //{
        //    std::cout << "U old in AB2 update " << i << '\n' << soln.u_old << '\n';
        //    std::cout << "U in AB2 update " << i << '\n' << soln.u << '\n';
        //}
        rhs.d_AB2_update(params, mesh, soln, op, boundaries);

        soln.u_old = soln.u;
        
        soln.u.head(soln.u.size() - 1) = thomasAlgorithm(matrix.CN_Matrix, rhs.d_AB2.head(rhs.d_AB2.size() - 1));
        X2 = thomasAlgorithm(matrix.CN_Matrix, rhs.q);
        
        soln.u(soln.u.size() - 1) = (rhs.d_AB2(rhs.d_AB2.size() - 1) - params.P * soln.u(0) - params.P * soln.u(soln.u.size() - 2)) /
            (params.Q + params.P * X2(0) + params.P * X2(X2.size() - 1));

        soln.u.head(soln.u.size() - 1) = soln.u.head(soln.u.size() - 1) + X2 * soln.u(soln.u.size() - 1);

        //std::cout << "U NEW in AB2 update " << i << '\n' << soln.u << '\n';

        writeToFile(file, soln.u);
    }

    file.close();
    
    std::cout << "Solution is "<< '\n' << soln.u << '\n';

    return 0;
}