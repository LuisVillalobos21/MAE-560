#include "Mesh.hpp"

Mesh::Mesh(double start, double end, double step)
{
    domain_start = start;
    domain_end = end;
    dx = step;
    number_points = static_cast<int>((domain_end + dx - domain_start) / dx);
    x = Eigen::VectorXd::LinSpaced(number_points, domain_start, domain_end);
}
