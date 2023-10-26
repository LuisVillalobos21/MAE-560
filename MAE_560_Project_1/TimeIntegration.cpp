#include "TimeIntegration.hpp"

void TimeIntegration::explicitEuler(const Parameters& params, const Mesh& mesh, SolutionVector& soln,
    const BoundaryConditions& boundaries, const Discretization& discrete)
{
    soln.u_new = soln.u + ((params.alpha * params.dt) / (mesh.dx * mesh.dx)) * (discrete.laplacian_matrix * soln.u + boundaries.u_bc);
    soln.u = soln.u_new;
}