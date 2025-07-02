#include <iostream>
#include <Eigen/Dense>
#include <chrono>
// #include <casadi/casadi.hpp>

using casadi_int = long long int;

// External C function declarations
extern "C" {
    int obj_fun(const double** arg, double** res, casadi_int* iw, double* w, int mem);
    int grad_fun(const double** arg, double** res, casadi_int* iw, double* w, int mem);
    int proj_fun(const double** arg, double** res, casadi_int* iw, double* w, int mem);

    int obj_fun_work(casadi_int*, casadi_int*, casadi_int*, casadi_int*);
    int grad_fun_work(casadi_int*, casadi_int*, casadi_int*, casadi_int*);
    int proj_fun_work(casadi_int*, casadi_int*, casadi_int*, casadi_int*);
}

/**
 * @brief CasADi PGD API class using static libraries
 * 
 * This class implements the Projected Gradient Descent (PGD) algorithm
 * using CasADi-generated static libraries for objective function,
 * gradient function, and projection function evaluations.
 * If libpgd_fun.a doesn't exist, it will automatically generate and build it.
 */
class PGD_API_s {
public:
    /**
     * @brief Constructor initializes the PGD solver with default parameters
     */
    PGD_API_s() {
        C[0] = 0.0; C[1] = 1.2;
        radius = 0.5;
    }

    /**
     * @brief Destructor
     */
    ~PGD_API_s() = default;

    /**
     * @brief Solve PGD optimization problem with Eigen interface
     * @param X_init Initial guess as Eigen::Vector2d
     */
    void solve(const Eigen::Vector2d& X_init) {
        double x_init_arr[2] = { X_init(0), X_init(1) };
        solve_internal(x_init_arr);
    }

    /**
     * @brief Display optimization results and timing information
     */
    void display_info() const {
        std::cout << "[PGD_API_s] total time: " << duration.count() << "s" << std::endl;
        std::cout << "[PGD_API_s] iter: " << iter << ", objective: " << FX
                  << ", X: (" << X[0] << ", " << X[1] << ")" << std::endl;
    }

    std::chrono::duration<double> duration;

private:
    /**
     * @brief Internal solver implementation
     * @param X_init Initial guess array
     */
    void solve_internal(const double X_init[2]) {
        auto tic = std::chrono::high_resolution_clock::now();

        // Initialize variables
        for (int i = 0; i < 2; i++) {
            X[i] = X_init[i];
            Y[i] = X[i];
        }

        evaluate_obj(Y, FX);
        FX_prev = FX;

        for (int i = 0; i < 2; i++) {
            Y_prev[i] = 0.0;
            grad_Y_prev[i] = 0.0;
        }

        tk = 1.0;
        qk = 1.0;
        ck = FX;

        eta = 0.4;
        del = 0.001;
        max_iter = 100;
        iter = 0;

        // Main optimization loop
        while (true) {
            iter++;

            evaluate_grad(Y, grad_Y);

            for (int i = 0; i < 2; i++) {
                sk[i] = Y[i] - Y_prev[i];
                rk[i] = grad_Y[i] - grad_Y_prev[i];
            }

            double numerator = sk[0]*rk[0] + sk[1]*rk[1];
            double denominator = rk[0]*rk[0] + rk[1]*rk[1];
            alpha_Y = std::abs(numerator / denominator);

            rho_Y = 0.5;
            back_iter = 0;

            // Backtracking line search for Y
            while (true) {
                back_iter++;
                for (int i = 0; i < 2; i++)
                    temp[i] = Y[i] - alpha_Y * grad_Y[i];
                evaluate_proj(temp, Z);

                alpha_Y *= rho_Y;
                evaluate_obj(Z, FZ);

                double diff = (ck - FZ) - del * squared_distance(Y, Z);
                if (diff >= 0 || back_iter > 10) break;
            }

            bool accepted_Z = (ck - FZ) >= del * squared_distance(Y, Z);
            if (accepted_Z) {
                for (int i = 0; i < 2; i++) X[i] = Z[i];
                FX = FZ;
            } else {
                evaluate_grad(X, grad_X);
                for (int i = 0; i < 2; i++) {
                    sk[i] = X[i] - Y_prev[i];
                    rk[i] = grad_X[i] - grad_Y_prev[i];
                }
                alpha_X = std::abs( (sk[0]*rk[0] + sk[1]*rk[1]) / (rk[0]*rk[0] + rk[1]*rk[1]) );

                rho_X = 0.5;
                mon_iter = 0;

                // Backtracking line search for X
                while (true) {
                    mon_iter++;
                    for (int i = 0; i < 2; i++)
                        temp[i] = X[i] - alpha_X * grad_X[i];
                    evaluate_proj(temp, V);

                    alpha_X *= rho_X;
                    evaluate_obj(V, FV);

                    double diff = (ck - FV) - del * squared_distance(Y, V);
                    if (diff >= 0 || mon_iter > 10) break;
                }

                if (FZ <= FV) {
                    for (int i = 0; i < 2; i++) X[i] = Z[i];
                    FX = FZ;
                } else {
                    for (int i = 0; i < 2; i++) X[i] = V[i];
                    FX = FV;
                }
            }

            tk_plus = (1 + sqrt(1 + 4 * tk * tk)) / 2.0;
            qk_plus = eta * qk + 1;
            ck_plus = (eta * qk * ck + FX) / qk_plus;

            if (pow(ck_plus - ck, 2) < 1e-6 || iter >= max_iter)
                break;

            for (int i = 0; i < 2; i++) {
                Y_prev[i] = Y[i];
                grad_Y_prev[i] = grad_Y[i];
                Y[i] = X[i] + tk / tk_plus * (Z[i] - X[i]) + (tk - 1) / tk_plus * (X[i] - X_prev[i]);
                X_prev[i] = X[i];
            }

            tk = tk_plus;
            qk = qk_plus;
            ck = ck_plus;
        }

        auto toc = std::chrono::high_resolution_clock::now();
        duration = toc - tic;
    }

    /**
     * @brief Evaluate objective function using CasADi
     * @param x Input point
     * @param result Output objective value
     */
    void evaluate_obj(const double x[2], double& result) {
        const double* arg[1] = {x};
        double* res[1] = {&result};
        casadi_int sz_arg, sz_res, sz_iw, sz_w;
        obj_fun_work(&sz_arg, &sz_res, &sz_iw, &sz_w);
        casadi_int iw[sz_iw];
        double w[sz_w];
        obj_fun(arg, res, iw, w, 0);
    }

    /**
     * @brief Evaluate gradient function using CasADi
     * @param x Input point
     * @param grad_out Output gradient
     */
    void evaluate_grad(const double x[2], double grad_out[2]) {
        const double* arg[1] = {x};
        double* res[1] = {grad_out};
        casadi_int sz_arg, sz_res, sz_iw, sz_w;
        grad_fun_work(&sz_arg, &sz_res, &sz_iw, &sz_w);
        casadi_int iw[sz_iw];
        double w[sz_w];
        grad_fun(arg, res, iw, w, 0);
    }

    /**
     * @brief Evaluate projection function using CasADi
     * @param input Input point
     * @param proj_out Output projected point
     */
    void evaluate_proj(const double input[2], double proj_out[2]) {
        const double* arg[3] = {input, C, &radius};
        double* res[1] = {proj_out};
        casadi_int sz_arg, sz_res, sz_iw, sz_w;
        // proj_fun_work_ptr(&sz_arg, &sz_res, &sz_iw, &sz_w);
        proj_fun_work(&sz_arg, &sz_res, &sz_iw, &sz_w);
        casadi_int iw[sz_iw];
        double w[sz_w];
        // proj_fun_ptr(arg, res, iw, w, 0);
        proj_fun(arg, res, iw, w, 0);
    }

    /**
     * @brief Calculate squared distance between two points
     * @param a First point
     * @param b Second point
     * @return Squared distance
     */
    double squared_distance(const double a[2], const double b[2]) const {
        return (a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]);
    }

    // Problem parameters
    double C[2], radius;
    
    // Optimization variables
    double X[2], Y[2], X_prev[2], Y_prev[2];
    double grad_Y[2], grad_Y_prev[2], grad_X[2];
    double Z[2], V[2], temp[2], sk[2], rk[2];
    double FX, FX_prev, FZ, FV;
    double tk, tk_plus, qk, qk_plus, ck, ck_plus;
    double eta, del;
    double alpha_Y, alpha_X;
    int rho_Y, rho_X;
    int iter, max_iter;
    int back_iter, mon_iter;
}; 