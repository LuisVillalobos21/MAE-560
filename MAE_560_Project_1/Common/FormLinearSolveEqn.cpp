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
	b = ((params.alpha * params.dt) / (2 * mesh.dx * mesh.dx)) * (discrete.laplacian_matrix * soln.u) +
		((params.alpha * params.dt) / (2 * mesh.dx * mesh.dx)) * (boundaries.u_bc + boundaries.u_bc) + soln.u;

	//std::cout << "The right hand side in function is" << b << '\n';
	//std::cout << "alpha is " << params.alpha << '\n';
	//std::cout << "dt is " << params.dt << '\n';
	//std::cout << "dx is " << mesh.dx << '\n';
	//std::cout << "laplacian is " << discrete.laplacian_matrix << '\n';
	//std::cout << "u is " << soln.u << '\n';
	//std::cout << "ubc vector is " << boundaries.u_bc << '\n';

}