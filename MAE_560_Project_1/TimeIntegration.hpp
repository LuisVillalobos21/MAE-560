#pragma once
#include <Eigen/Dense>
#include "Parameters.hpp"
#include "Mesh.hpp"
#include "SolutionVector.hpp"
#include "BoundaryConditions.hpp"
#include "Discretization.hpp"

struct TimeIntegration
{
    void explicitEuler(const Parameters& params, const Mesh& mesh, SolutionVector& soln,
        const BoundaryConditions& boundaries, const Discretization& discrete);
};