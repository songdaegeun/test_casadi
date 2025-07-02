#include <casadi/casadi.hpp>
#include <iostream>

using namespace casadi;
using namespace std;

int main() {

    MX x = MX::sym("x", 2);  // [x1, x2]

    // 목적 함수: 15*(x1^2-1)^2 + 1*(x2^2-2)^2 + 4*x1*x2 + x1 + x2 
    MX f = 15.0*pow(x(0)*x(0)-1, 2) + 1.0*pow(x(1)*x(1) - 2, 2) + 4.0*x(0)*x(1) + x(0) + x(1);


    Function f_eval = Function("f_eval", {x}, {f});
    std::vector<double> x_val_1 = {0.499946, 1.19269};
    DM result_1 = f_eval(DM(x_val_1))[0];
    std::cout << "[Ipopt]: f(" << x_val_1[0] << ", " << x_val_1[1] << ") = " << result_1 << std::endl;


    
    std::vector<double> x_val_2 = {0.499947, 1.19269}; // 1e-6 -> 8.48031
    // std::vector<double> x_val_2 = {0.499747, 1.1841}; // 1e-0 -> 8.47034
    DM result_2 = f_eval(DM(x_val_2))[0];
    std::cout << "[PGD]: f(" << x_val_2[0] << ", " << x_val_2[1] << ") = " << result_2 << std::endl;
    
    // 제약식: (x1)^2 + (x2-1.2)^2 <= 0.5^2
    MX g1 = pow(x(0), 2) + pow(x(1) - 1.2, 2);

    Function g_eval = Function("g_eval", {x}, {g1});
    DM result_3 = g_eval(DM(x_val_1))[0];
    std::cout << "[Ipopt]: g(" << x_val_1[0] << ", " << x_val_1[1] << ") = " << result_3 << std::endl;

    DM result_4 = g_eval(DM(x_val_2))[0];
    std::cout << "[PGD]: g(" << x_val_2[0] << ", " << x_val_2[1] << ") = " << result_4 << std::endl;


    std::vector<MX> g = {g1};

    // cout << "g: " << g << endl;
    // cout << "vertcat(g): " << vertcat(g) << endl;

    MXDict nlp;
    nlp["x"] = x;
    nlp["f"] = f;
    nlp["g"] = vertcat(g);

    Dict opts;
    // opts["ipopt.max_iter"] = 10000; 
    // opts["ipopt.tol"] = 1e-4;
    // opts["ipopt.tol"] = 1e-8;
    // opts["ipopt.tol"] = 1e-14;
    opts["ipopt.tol"] = 1e-4;
    opts["ipopt.dual_inf_tol"] = 1e-4;
    opts["ipopt.constr_viol_tol"] = 1e-4;
    opts["ipopt.compl_inf_tol"] = 1e-4;
    // opts["ipopt.linear_solver"] = "ma57";
// ==============

    auto tic1 = std::chrono::high_resolution_clock::now();
    Function solver = nlpsol("solver", "ipopt", nlp, opts);
    auto tac1 = std::chrono::high_resolution_clock::now();

    std::vector<double> x0 = {0.2, 1};
    // std::vector<double> x0 = {-0.07, -1.45};

    // 변수 경계 조건
    std::vector<double> lbx = {-inf, -inf};
    std::vector<double> ubx = {inf, inf};

    // 제약 조건 경계
    std::vector<double> lbg = {-inf};
    std::vector<double> ubg = {0.5*0.5};

    DMDict arg;
    arg["x0"] = DM(x0);
    arg["lbx"] = DM(lbx);
    arg["ubx"] = DM(ubx);
    arg["lbg"] = DM(lbg);
    arg["ubg"] = DM(ubg);

    auto tic2 = std::chrono::high_resolution_clock::now();
    DMDict res = solver(arg);
    auto tac2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration1 = tac1 - tic1;
    std::chrono::duration<double> duration2 = tac2 - tic2;
    std::cout << "Setting Time taken: " << duration1.count() << " seconds" << std::endl;
    std::cout << "Solving Time taken: " << duration2.count() << " seconds" << std::endl;

    // 결과 출력
    cout << "optimal X: " << res.at("x") << endl;
    cout << "optimal F: " << res.at("f") << endl;

    return 0;
}