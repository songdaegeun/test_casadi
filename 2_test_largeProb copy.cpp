#include <casadi/casadi.hpp>
using namespace casadi;
int main() {
    int N = 500;
    MX x = MX::sym("x", N);

    MX f = 0;
    for (int i = 0; i < N; ++i)
        f += pow(x(i) - i, 4) + sin(x(i)) + exp(0.1 * x(i));

    std::vector<MX> g;
    for (int i = 0; i < N-1; ++i)
        g.push_back(pow(x(i), 2) + pow(x(i+1), 2) - 1.0);

    MXDict nlp;
    nlp["x"] = x;
    nlp["f"] = f;
    nlp["g"] = vertcat(g);

    Dict opts;
    opts["ipopt.max_iter"] = 10000; 
    opts["ipopt.tol"] = 1e-10;

    Function solver = nlpsol("solver", "ipopt", nlp, opts);

    std::vector<double> x0(N, 0.0);
    std::vector<double> lbg(N-1, 0.0), ubg(N-1, 0.0);

    DMDict arg;
    arg["x0"] = DM(x0);
    arg["lbg"] = DM(lbg);
    arg["ubg"] = DM(ubg);

    DMDict res = solver(arg);
    std::cout << "최적해: " << res.at("x") << std::endl;
    std::cout << "최적값: " << res.at("f") << std::endl;
    return 0;
}