// Rosenbrock problem: https://web.casadi.org/blog/opti/

#include "rosenbrock.h"

int main() {
    casadi::Opti opti;
    
    casadi::MX x = opti.variable();
    casadi::MX y = opti.variable();
    casadi::MX f = pow(1 - x, 2) + pow(y - x * x, 2);

    opti.minimize(f);
    casadi::DM X = casadi::DM::linspace(0, 1.5, 100);
    casadi::DM Y = casadi::DM::linspace(-0.5, 1.5, 100);
    auto [XX, YY] = meshgrid(X, Y);
    casadi::DM ZZ = (1-XX)*(1-XX) + (YY-XX*XX)*(YY-XX*XX);
    save_contour(XX, YY, ZZ);

    // casadi::MX r = opti.parameter();
    // opti.subject_to(x*x+y*y<=r);
    // save_constraint_1(casadi::DM(r).scalar());

    opti.subject_to(x*x+y*y<=1);
    // save_constraint_1();
    opti.subject_to(y>=x);
    save_constraint_2();
    
    opti.solver("ipopt");

    // std::vector<double> r_values;
    // std::vector<double> f_values;
    // casadi::DM r_list = casadi::DM::linspace(1.0, 3.0, 25);
    // for (int i = 0; i < r_list.size1(); ++i) {
    //   double r_val = static_cast<double>(r_list(i));
    //   opti.set_value(r, r_val);
    //   casadi::OptiSol sol = opti.solve();
    //   double f_val = static_cast<double>(sol.value(f));
    //   r_values.push_back(r_val);
    //   f_values.push_back(f_val);
    // }

    casadi::OptiSol sol = opti.solve();

    double x_opt = sol.value(x).scalar();
    double y_opt = sol.value(y).scalar();
    save_optimal_solution(x_opt, y_opt);

    return 0;
}

std::pair<casadi::DM, casadi::DM> meshgrid(const casadi::DM& x, const casadi::DM& y) {
  int nx = x.size1();
  int ny = y.size1();
  casadi::DM X(ny, nx);
  casadi::DM Y(ny, nx);

  for (int i = 0; i < ny; ++i) {
    for (int j = 0; j < nx; ++j) {
      X(i, j) = x(j);
      Y(i, j) = y(i);
    }
  }
  return {X, Y};
}

void save_dm(const casadi::DM& data, const char* name, mat_t* matfp) {
  size_t dims[2] = {(size_t)data.size1(), (size_t)data.size2()};
  matvar_t* matvar = Mat_VarCreate(name, MAT_C_DOUBLE, MAT_T_DOUBLE,
                                    2, dims, (void*)data.ptr(), 0);
  if (!matvar) throw std::runtime_error("Failed to create matvar");
  Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_NONE);
  Mat_VarFree(matvar);
}

void save_contour(const casadi::DM& XX, const casadi::DM& YY, const casadi::DM& ZZ) {
  mat_t* matfp = Mat_CreateVer("../mfiles/mat/rosenbrock_contour.mat", NULL, MAT_FT_MAT5);
  if (!matfp) throw std::runtime_error("Failed to create .mat file");

  save_dm(XX, "XX", matfp);
  save_dm(YY, "YY", matfp);
  save_dm(ZZ, "ZZ", matfp);

  Mat_Close(matfp);
}

void save_constraint_1(double r) {
  mat_t* matfp = Mat_CreateVer("../mfiles/mat/rosenbrock_constraint_1.mat", NULL, MAT_FT_MAT5);
  if (!matfp) throw std::runtime_error("Failed to create .mat file");

  int N = 500;
  std::vector<double> circle_x(N), circle_y(N);
  for (int i = 0; i < N; ++i) {
      double theta = 2 * M_PI * i / (N - 1);
      circle_x[i] = r * std::cos(theta);
      circle_y[i] = r * std::sin(theta);
  }

  size_t dims1[2] = {1, (size_t)N};
  matvar_t* var_cx = Mat_VarCreate("circle_x", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims1, circle_x.data(), 0);
  matvar_t* var_cy = Mat_VarCreate("circle_y", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims1, circle_y.data(), 0);
  if (!var_cx || !var_cy) throw std::runtime_error("Failed to create matvar");
  Mat_VarWrite(matfp, var_cx, MAT_COMPRESSION_NONE);
  Mat_VarWrite(matfp, var_cy, MAT_COMPRESSION_NONE);
  Mat_VarFree(var_cx);
  Mat_VarFree(var_cy);

  Mat_Close(matfp);
}

void save_constraint_2() {
  mat_t* matfp = Mat_CreateVer("../mfiles/mat/rosenbrock_constraint_2.mat", NULL, MAT_FT_MAT5);
  if (!matfp) throw std::runtime_error("Failed to create .mat file");

  int M = 200;
  std::vector<double> line_x(M), line_y(M);
  for (int i = 0; i < M; ++i) {
      line_x[i] = 0.0 + 1.5 * i / (M - 1); // 0 ~ 1.5
      line_y[i] = line_x[i];
  }

  size_t dims2[2] = {1, (size_t)M};
  matvar_t* var_lx = Mat_VarCreate("line_x", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims2, line_x.data(), 0);
  matvar_t* var_ly = Mat_VarCreate("line_y", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims2, line_y.data(), 0);
  if (!var_lx || !var_ly) throw std::runtime_error("Failed to create matvar");
  Mat_VarWrite(matfp, var_lx, MAT_COMPRESSION_NONE);
  Mat_VarWrite(matfp, var_ly, MAT_COMPRESSION_NONE);
  Mat_VarFree(var_lx);
  Mat_VarFree(var_ly);

  Mat_Close(matfp);
}

void save_optimal_solution( double x_opt, double y_opt) {
  mat_t* matfp = Mat_CreateVer("../mfiles/mat/rosenbrock_optimal_solution.mat", NULL, MAT_FT_MAT5);
  if (!matfp) throw std::runtime_error("Failed to create .mat file");

  casadi::DM xopt_dm = casadi::DM(std::vector<std::vector<double>>{{x_opt}});
  casadi::DM yopt_dm = casadi::DM(std::vector<std::vector<double>>{{y_opt}});
  save_dm(xopt_dm, "x_opt", matfp);
  save_dm(yopt_dm, "y_opt", matfp);

  Mat_Close(matfp);
}
