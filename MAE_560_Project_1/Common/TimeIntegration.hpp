#pragma once
#include <Eigen/Dense>
#include "Parameters.hpp"
#include "Mesh.hpp"
#include "SolutionVector.hpp"
#include "BoundaryConditions.hpp"
#include "Discretization.hpp"
#include "FormLinearSolveEqn.hpp"

struct TimeIntegration
{
    void explicitEuler(const Parameters& params, const Mesh& mesh, SolutionVector& soln,
        const BoundaryConditions& boundaries, const Discretization& discrete);
};

struct CNTimeIntegration
{
    void thomasAlgorithm(const CNDiscretization& discrete, const FormLinearSolveEqn& rhs, SolutionVector& soln);
};