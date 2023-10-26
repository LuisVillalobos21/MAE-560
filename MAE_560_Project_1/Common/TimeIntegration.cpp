#include "TimeIntegration.hpp"

void TimeIntegration::explicitEuler(const Parameters& params, const Mesh& mesh, SolutionVector& soln,
    const BoundaryConditions& boundaries, const Discretization& discrete)
{
    soln.u_new = soln.u + ((params.alpha * params.dt) / (mesh.dx * mesh.dx)) * (discrete.laplacian_matrix * soln.u + boundaries.u_bc);
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