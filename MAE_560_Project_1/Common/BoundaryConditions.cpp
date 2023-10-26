#include "BoundaryConditions.hpp"

BoundaryConditions::BoundaryConditions(int number_points, const Parameters& params, const Mesh& mesh,
    const SolutionVector& soln)
{
    u_bc = Eigen::VectorXd::Zero(number_points);
    bc_rhs = Eigen::VectorXd::Zero(2);
    robin_matrix = createRobinLeftBCMatrix(params, mesh);
    calcRobinValue(params, mesh, soln);
    u_bc(u_bc.size() - 1) = u_bc_L;
}

Eigen::MatrixXd BoundaryConditions::createRobinLeftBCMatrix(const Parameters& params, const Mesh& mesh)
{
    Eigen::MatrixXd bc_matrix = Eigen::MatrixXd::Zero(2, 2);
    bc_matrix(0, 0) = 2;
    bc_matrix(0, 1) = -1;
    bc_matrix(1, 0) = -(2 * params.a * mesh.dx) / params.b;
    bc_matrix(1, 1) = 1;
    return bc_matrix;
}

void BoundaryConditions::calcRobinValue(const Parameters& params, const Mesh& mesh, const SolutionVector& soln)
{
    bc_rhs(0) = soln.u(0);
    bc_rhs(1) = soln.u(0) - (2 * params.c * mesh.dx) / params.b;

    Eigen::VectorXd f_bc = robin_matrix.partialPivLu().solve(bc_rhs);
    u_bc(0) = f_bc(0);
}