#include "FormLinearSolveEqn.hpp"

FormLinearSolveEqn::FormLinearSolveEqn(const int number_points, const Parameters& params, const Mesh& mesh,
	const SolutionVector& soln, const BoundaryConditions& boundaries, const CNDiscretization& discrete)
{
	b = Eigen::VectorXd::Zero(number_points);
	calcLinearSolveRHS(params, mesh, soln, boundaries, discrete);
}

void FormLinearSolveEqn::calcLinearSolveRHS(const Parameters& params, const Mesh& mesh,
	const SolutionVector& soln, const BoundaryConditions& boundaries, const CNDiscretization& discrete)
{
	b = ((params.alpha * params.dt) / (2 * mesh.dx * mesh.dx)) * (discrete.laplacian_matrix * soln.u) + boundaries.u_bc + soln.u;
}