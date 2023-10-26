#pragma once
#include <Eigen/Dense>
#include "Parameters.hpp"
#include "Mesh.hpp"
#include "SolutionVector.hpp"
#include "BoundaryConditions.hpp"
#include "Discretization.hpp"
#include <iostream>

struct FormLinearSolveEqn
{
    Eigen::VectorXd b;

    FormLinearSolveEqn(const int number_points, const Parameters& params, const Mesh& mesh,
        const SolutionVector& soln, const BoundaryConditions& boundaries, const CNDiscretization& discrete);

    void calcLinearSolveRHS(const Parameters& params, const Mesh& mesh, const SolutionVector& soln, 
        const BoundaryConditions& boundaries, const CNDiscretization& discrete);

};
