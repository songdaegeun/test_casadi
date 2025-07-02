#ifndef ROSENBROCK_H
#define ROSENBROCK_H

#include <casadi/casadi.hpp>
#include <matio.h>

std::pair<casadi::DM, casadi::DM> meshgrid(const casadi::DM& x, const casadi::DM& y);
void save_dm(const casadi::DM& data, const char* name, mat_t* matfp);
void save_contour(const casadi::DM& XX, const casadi::DM& YY, const casadi::DM& ZZ);
void save_optimal_solution( double x_opt, double y_opt);
void save_constraint_1(double r = 1);
void save_constraint_2();

#endif