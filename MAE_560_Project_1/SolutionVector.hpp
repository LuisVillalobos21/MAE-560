#pragma once
#include <Eigen/Dense>

struct SolutionVector
{
    Eigen::VectorXd u;
    Eigen::VectorXd u_new;
    SolutionVector(int number_points);
};