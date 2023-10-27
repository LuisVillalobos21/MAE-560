#pragma once
#include <Eigen/Dense>

struct SolutionVector
{
    Eigen::VectorXd u;
    Eigen::VectorXd u_new;
    Eigen::VectorXd u_old;
    SolutionVector(int number_points);
};