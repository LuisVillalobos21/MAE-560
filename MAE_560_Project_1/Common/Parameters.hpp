#pragma once
#include <Eigen/Dense>

struct Parameters
{
    double alpha;
    double dt;
    double a;
    double b;
    double c;
};

struct WaveParameters
{
    double c;
    double dt;
};