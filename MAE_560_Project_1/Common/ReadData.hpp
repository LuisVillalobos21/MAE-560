#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <Eigen/Dense>

Eigen::VectorXd readCSVtoEigenVector(const std::string& filename);