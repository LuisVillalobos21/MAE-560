#include "BoundaryConditions.hpp"

// ROBIN BC AND DIRECHLET FOR HEAT EQN
BoundaryConditions::BoundaryConditions(int number_points, const Parameters& params, const Mesh& mesh,
    const SolutionVector& soln)
{
    u_bc = Eigen::VectorXd::Zero(number_points);
    calcRobinValue(params, mesh, soln);
    u_bc(u_bc.size() - 1) = ((params.alpha * params.dt) / (mesh.dx * mesh.dx)) * u_bc_R;

    u_bc_exp = Eigen::VectorXd::Zero(number_points);
    calcRobinValueExp(params, mesh, soln);
    u_bc_exp(u_bc_exp.size() - 1) = ((params.alpha * params.dt) / ( mesh.dx * mesh.dx)) * u_bc_R;
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
    u_bc(0) = ((params.alpha * params.dt) / (2 * mesh.dx * mesh.dx)) * ((params.c - (params.b / mesh.dx) * soln.u(0)) / (params.a - (params.b / mesh.dx)) + params.c / (params.a - (params.b / mesh.dx)));
}

void BoundaryConditions::calcRobinValueExp(const Parameters& params, const Mesh& mesh, const SolutionVector& soln)
{
    u_bc_exp(0) = ((params.alpha * params.dt) / (mesh.dx * mesh.dx)) * ((params.c - (params.b / mesh.dx) * soln.u(0)) / (params.a - (params.b / mesh.dx)));
}

// PERIODIC BC FOR WAVE EQN LAX WENDROFF
PeriodicBoundaryConditions::PeriodicBoundaryConditions(const WaveParameters& params, const PeriodicMesh& mesh, const SolutionVector& soln)
{
    u_bc_diff = Eigen::VectorXd::Zero(mesh.number_points);
    u_bc_lap = Eigen::VectorXd::Zero(mesh.number_points);
    calcPeriodicValue(soln,mesh);
}

void PeriodicBoundaryConditions::calcPeriodicValue(const SolutionVector& soln, const PeriodicMesh& mesh)
{
    //u_bc_diff(0) = -soln.u(soln.u.size() - 1);    
    //u_bc_diff(u_bc_diff.size() - 1) = soln.u(0);

    //u_bc_lap(0) = soln.u(soln.u.size() - 1);
    //u_bc_lap(u_bc_lap.size() - 1) = soln.u(0);

    u_bc_diff(0) = 0;
    u_bc_diff(u_bc_diff.size() - 1) = 0;

    u_bc_lap(0) = 0;
    u_bc_lap(u_bc_lap.size() - 1) = 0;
}

// PERIODIC BC FOR WAVE EQN AB2
AB2PeriodicBoundaryConditions::AB2PeriodicBoundaryConditions(const WaveParameters& params, const PeriodicMesh& mesh, const SolutionVector& soln)
{
    u_bc_diff1 = Eigen::VectorXd::Zero(mesh.number_points);
    u_bc_diff2 = Eigen::VectorXd::Zero(mesh.number_points);
    calcPeriodicValueAB2(soln, mesh);
}

void AB2PeriodicBoundaryConditions::calcPeriodicValueAB2(const SolutionVector& soln, const PeriodicMesh& mesh)
{
    u_bc_diff1(0) = -soln.u(soln.u.size() - 1);
    u_bc_diff1(u_bc_diff1.size() - 1) = soln.u(0);

    u_bc_diff2(0) = -soln.u_old(soln.u.size() - 1);
    u_bc_diff2(u_bc_diff2.size() - 1) = soln.u_old(0);
}