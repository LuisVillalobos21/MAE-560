#include "SolutionVector.hpp"

SolutionVector::SolutionVector(int number_points)
{
    u = Eigen::VectorXd::Zero(number_points);
    u_new = Eigen::VectorXd::Zero(number_points);
}