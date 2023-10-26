#pragma once
#include <Eigen/Dense>

struct Mesh
{
    double dx;
    double domain_start;
    double domain_end;
    int number_points;
    Eigen::VectorXd x;

    Mesh(double start, double end, double step);
};
