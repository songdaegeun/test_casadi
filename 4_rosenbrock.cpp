// Rosenbrock problem: https://web.casadi.org/blog/opti/
// opti: https://github.com/casadi/casadi/blob/develop/docs/examples/cplusplus/race_car.cpp

#include <casadi/casadi.hpp>
#include <vector>

std::pair<casadi::DM, casadi::DM> meshgrid(const casadi::DM& x, const casadi::DM& y);

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
    casadi::MX ZZ = (1-XX)*(1-XX) + (YY-XX*XX)*(YY-XX*XX);
    casadi::MX Z = casadi::MX::reshape(ZZ, 100, 100);

    // // Create Matlab script to plot the solution
    // std::ofstream file;
    // std::string filename = "rosenbrock_result.m";
    // file.open(filename.c_str());
    // file << "% Results file from " __FILE__ << std::endl;
    // file << "% Generated " __DATE__ " at " __TIME__ << std::endl;
    // file << std::endl;
    
    // // Save results to file
    // file << "x_opt = " << x_opt << ";" << std::endl;
    // file << "y_opt = " << y_opt << ";" << std::endl;
    // file << "figure;\n";
    // file << "plot(x_opt, y_opt, 'ro', 'MarkerSize', 10, 'LineWidth', 2);\n";
    // file << "xlabel('x');\n";
    // file << "ylabel('y');\n";
    // file << "title('Rosenbrock solution');\n";
    // file << "grid on;\n";
    // file << "axis equal;\n";
    // file.close();

    // std::cout << "Rosenbrock solution saved to " << filename << std::endl;
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