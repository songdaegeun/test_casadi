// Rosenbrock problem: https://web.casadi.org/blog/opti/
// opti: https://github.com/casadi/casadi/blob/develop/docs/examples/cplusplus/race_car.cpp

#include <casadi/casadi.hpp>
#include <matio.h>

std::pair<casadi::DM, casadi::DM> meshgrid(const casadi::DM& x, const casadi::DM& y);
void save_matlab_data(const casadi::DM& XX, const casadi::DM& YY, const casadi::DM& ZZ, double x_opt, double y_opt);

int main() {
    casadi::Opti opti;
    
    casadi::MX x = opti.variable();
    casadi::MX y = opti.variable();
    // 목적 함수: Rosenbrock
    opti.minimize(pow(1 - x, 2) + pow(y - x * x, 2));

    opti.solver("ipopt");

    casadi::OptiSol sol = opti.solve();

    double x_opt = sol.value(x).scalar();
    double y_opt = sol.value(y).scalar();

    // std::cout << "x_opt: " << x_opt << std::endl;
    // std::cout << "y_opt: " << y_opt << std::endl;

    casadi::DM X = casadi::DM::linspace(-2, 2, 100);
    casadi::DM Y = casadi::DM::linspace(-2, 2, 100);
    auto [XX, YY] = meshgrid(X, Y);
    casadi::DM ZZ = (1-XX)*(1-XX) + (YY-XX*XX)*(YY-XX*XX);

    // std::cout << "XX.shape: " << XX.size1() << " " << XX.size2() << std::endl;
    // std::cout << "YY.shape: " << YY.size1() << " " << YY.size2() << std::endl;
    // std::cout << "ZZ.shape: " << ZZ.size1() << " " << ZZ.size2() << std::endl;

    save_matlab_data(XX, YY, ZZ, x_opt, y_opt);

    // // 파일 열기
    // std::ofstream file("../mfiles/rosenbrock_result.m");
    // file << "% Rosenbrock contour plot\n";

    // // 4. 나머지 시각화 코드
    // file << "figure;\n";
    // file << "contour(XX, YY, ZZ, 100);\n";
    // file << "colormap('viridis');\n";
    // file << "colorbar;\n";
    // file << "hold on;\n";
    // file << "plot(" << x_opt << ", " << y_opt << ", 'ro', 'MarkerSize', 10, 'LineWidth', 2);\n";
    // file << "xlabel('x'); ylabel('y'); title('Rosenbrock Contour');\n";
    // file << "grid on;\n";
    // file.close();

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