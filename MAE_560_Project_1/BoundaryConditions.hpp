#pragma once
#include <Eigen/Dense>
#include "Parameters.hpp"
#include "Mesh.hpp"
#include "SolutionVector.hpp"

struct BoundaryConditions
{
    Eigen::VectorXd u_bc;
    Eigen::MatrixXd robin_matrix;
    Eigen::VectorXd bc_rhs;
    double u_bc_L = 1;

    BoundaryConditions(int number_points, const Parameters& params, const Mesh& mesh,
        const SolutionVector& soln);

    Eigen::MatrixXd createRobinLeftBCMatrix(const Parameters& params, const Mesh& mesh);

    void calcRobinValue(const Parameters& params, const Mesh& mesh, const SolutionVector& soln);
};