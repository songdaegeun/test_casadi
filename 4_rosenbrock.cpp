// Rosenbrock problem: https://web.casadi.org/blog/opti/
// opti: https://github.com/casadi/casadi/blob/develop/docs/examples/cplusplus/race_car.cpp

#include <casadi/casadi.hpp>
#include <matio.h>

std::pair<casadi::DM, casadi::DM> meshgrid(const casadi::DM& x, const casadi::DM& y);
void save_matlab_data(const casadi::DM& XX, const casadi::DM& YY, const casadi::DM& ZZ, double x_opt, double y_opt);
void make_matlab_script( );

int main() {
    casadi::Opti opti;
    
    casadi::MX x = opti.variable();
    casadi::MX y = opti.variable();
    // 목적 함수: Rosenbrock
    opti.minimize(pow(1 - x, 2) + pow(y - x * x, 2));
    // opti.subject_to(x*x+y*y<=1);

    opti.solver("ipopt");

    casadi::OptiSol sol = opti.solve();

    double x_opt = sol.value(x).scalar();
    double y_opt = sol.value(y).scalar();

    casadi::DM X = casadi::DM::linspace(-2, 2, 100);
    casadi::DM Y = casadi::DM::linspace(-2, 2, 100);
    auto [XX, YY] = meshgrid(X, Y);
    casadi::DM ZZ = (1-XX)*(1-XX) + (YY-XX*XX)*(YY-XX*XX);

    save_matlab_data(XX, YY, ZZ, x_opt, y_opt);
    // make_matlab_script();

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

void save_matlab_data(const casadi::DM& XX, const casadi::DM& YY, const casadi::DM& ZZ, double x_opt, double y_opt) {
  mat_t* matfp = Mat_CreateVer("../mfiles/rosenbrock_data.mat", NULL, MAT_FT_MAT5);
  if (!matfp) throw std::runtime_error("Failed to create .mat file");

  auto save_dm = [&](const casadi::DM& data, const char* name) {
    size_t dims[2] = {(size_t)data.size1(), (size_t)data.size2()};
    matvar_t* matvar = Mat_VarCreate(name, MAT_C_DOUBLE, MAT_T_DOUBLE,
                                     2, dims, (void*)data.ptr(), 0);
    if (!matvar) throw std::runtime_error("Failed to create matvar");
    Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_NONE);
    Mat_VarFree(matvar);
  };

  save_dm(XX, "XX");
  save_dm(YY, "YY");
  save_dm(ZZ, "ZZ");

  casadi::DM xopt_dm = casadi::DM(std::vector<std::vector<double>>{{x_opt}});
  casadi::DM yopt_dm = casadi::DM(std::vector<std::vector<double>>{{y_opt}});
  save_dm(xopt_dm, "x_opt");
  save_dm(yopt_dm, "y_opt");

  Mat_Close(matfp);
}

void make_matlab_script( ) {
  std::ofstream file("../mfiles/rosenbrock_result.m");
  file << "% Rosenbrock contour plot" << std::endl;
  file << "% Load precomputed data from CasADi C++ output" << std::endl;
  file << "load('rosenbrock_data.mat');  % Loads XX, YY, ZZ, x_opt, y_opt" << std::endl;

  file << "figure;" << std::endl;
  file << "contour(XX, YY, ZZ, 100);" << std::endl;
  file << "colormap('parula');" << std::endl;
  file << "colorbar;" << std::endl;
  file << "hold on;" << std::endl;
  file << "plot(x_opt, y_opt, 'ro', 'MarkerSize', 10, 'LineWidth', 2);" << std::endl;
  file << "xlabel('x'); ylabel('y'); title('Rosenbrock Contour');" << std::endl;
  file << "grid on;" << std::endl;
  file.close();
}