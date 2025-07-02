#include <casadi/casadi.hpp>
#include <iostream>

using namespace casadi;
using namespace std;

int main() {
    MX x = MX::sym("x", 2);  // [x1, x2]

    // 목적 함수: (x1-1)^2 + (x2-2)^2
    MX f = pow(x(0) - 1, 2) + pow(x(1) - 2, 2);

    // 제약식: x1^2 + x2^2 <= 1, x1 + x2 <= 1
    MX g1 = pow(x(0), 2) + pow(x(1), 2);
    // MX g2 = x(0) + x(1);
    // std::vector<MX> g = {g1, g2};
    std::vector<MX> g = {g1};

    // Print g and vertcat(g)
    cout << "g: " << g << endl;
    cout << "vertcat(g): " << vertcat(g) << endl;

    MXDict nlp;
    nlp["x"] = x;
    nlp["f"] = f;
    nlp["g"] = vertcat(g);

    Function solver = nlpsol("solver", "ipopt", nlp);

    std::vector<double> x0 = {0.5, 0.5};

    std::vector<double> lbx = {0.0, 0.0};
    std::vector<double> ubx = {inf, inf};

    // std::vector<double> lbg = {-inf, -inf};
    // std::vector<double> ubg = {1.0, 1.0};
    std::vector<double> lbg = {-inf};
    std::vector<double> ubg = {1.0};

    DMDict arg;
    arg["x0"] = DM(x0);
    arg["lbx"] = DM(lbx);
    arg["ubx"] = DM(ubx);
    arg["lbg"] = DM(lbg);
    arg["ubg"] = DM(ubg);

    DMDict res = solver(arg);

    cout << "최적해: " << res.at("x") << endl;
    cout << "최적값: " << res.at("f") << endl;

    return 0;
}