#include "casadi_api_a_test.hpp"

int main () {
    Eigen::Vector2d X_init(0.2, 1.0); 

    PGD_API_s casadi_pgd_s;
    auto tic4 = std::chrono::high_resolution_clock::now();
    casadi_pgd_s.solve(X_init);
    auto toc4 = std::chrono::high_resolution_clock::now();
    casadi_pgd_s.duration = toc4 - tic4;
    casadi_pgd_s.display_info();  
}