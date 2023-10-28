#pragma once
#include <Eigen/Dense>
#include "Parameters.hpp"
#include "Mesh.hpp"
#include "SolutionVector.hpp"
#include <iostream>
struct BoundaryConditions
{
    Eigen::VectorXd u_bc;
    Eigen::MatrixXd robin_matrix;
    Eigen::VectorXd bc_rhs;

    Eigen::VectorXd u_bc_exp;
    double u_bc_R = 1;

    BoundaryConditions(int number_points, const Parameters& params, const Mesh& mesh,
        const SolutionVector& soln);

    Eigen::MatrixXd createRobinLeftBCMatrix(const Parameters& params, const Mesh& mesh);

    void calcRobinValue(const Parameters& params, const Mesh& mesh, const SolutionVector& soln);

    void calcRobinValueExp(const Parameters& params, const Mesh& mesh, const SolutionVector& soln);
};

struct PeriodicBoundaryConditions
{
    Eigen::VectorXd u_bc_diff;
    Eigen::VectorXd u_bc_lap;

    PeriodicBoundaryConditions(const WaveParameters& params, const PeriodicMesh& mesh,const SolutionVector& soln);

    void calcPeriodicValue(const SolutionVector& soln, const PeriodicMesh& mesh);

};


struct AB2PeriodicBoundaryConditions
{
    Eigen::VectorXd u_bc_diff1;
    Eigen::VectorXd u_bc_diff2;

    AB2PeriodicBoundaryConditions(const WaveParameters& params, const PeriodicMesh& mesh, const SolutionVector& soln);
    void calcPeriodicValueAB2(const SolutionVector& soln, const PeriodicMesh& mesh);
};
