#include "TimeIntegration.hpp"

void TimeIntegration::explicitEuler(const Parameters& params, const Mesh& mesh, SolutionVector& soln,
    const BoundaryConditions& boundaries, const Discretization& discrete)
{
    soln.u_new = soln.u + ((params.alpha * params.dt) / (mesh.dx * mesh.dx)) * (discrete.laplacian_matrix * soln.u) + boundaries.u_bc_exp;
    soln.u = soln.u_new;
}

void CNTimeIntegration::thomasAlgorithm(const CNDiscretization& discrete, const FormLinearSolveEqn& rhs, SolutionVector& soln)
{
    int n = soln.u.size();

    Eigen::VectorXd a(n), b(n), c(n), d = rhs.b;

    for (int i = 0; i < n; i++) {
        b(i) = discrete.cn_matrix(i, i);
        if (i > 0) a(i) = discrete.cn_matrix(i, i - 1);
        if (i < n - 1) c(i) = discrete.cn_matrix(i, i + 1);
    }

    for (int i = 1; i < n; i++) {
        double m = a(i) / b(i - 1);
        b(i) = b(i) - m * c(i - 1);
        d(i) = d(i) - m * d(i - 1);
    }

    soln.u(n - 1) = d(n - 1) / b(n - 1);
    for (int i = n - 2; i >= 0; i--) {
        soln.u(i) = (d(i) - c(i) * soln.u(i + 1)) / b(i);
    }
}

void LWTimeIntegration::laxWendroff(const WaveParameters& params, const PeriodicMesh& mesh, SolutionVector& soln,
    const PeriodicBoundaryConditions& boundaries, const LWDiscretization& discrete)
{
    soln.u_new = soln.u + ((- 0.5 * params.c * params.dt) / mesh.dx) * (discrete.differentiaion_matrix * soln.u + boundaries.u_bc_diff) +
        ((0.5 * params.c * params.c * params.dt * params.dt) / (mesh.dx * mesh.dx)) * (discrete.laplacian_matrix * soln.u + boundaries.u_bc_lap);
    soln.u = soln.u_new;
}

void AB2TimeIntegration::AB2(const WaveParameters& params, const PeriodicMesh& mesh, SolutionVector& soln,
    const AB2PeriodicBoundaryConditions& boundaries, const LWDiscretization& discrete)
{
    soln.u_new = soln.u + ((.75 * -params.c * params.dt) / mesh.dx) * (discrete.differentiaion_matrix * soln.u + boundaries.u_bc_diff1) +
        ((.25 * params.c * params.dt) / mesh.dx) * (discrete.differentiaion_matrix * soln.u_old + boundaries.u_bc_diff2);
    soln.u_old = soln.u;
    soln.u = soln.u_new;
}