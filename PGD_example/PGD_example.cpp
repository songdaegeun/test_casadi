#include <iostream>
#include <chrono>
#include <cmath>
#include </usr/include/eigen3/Eigen/Dense>
#include </usr/include/eigen3/unsupported/Eigen/KroneckerProduct> 
using Eigen::MatrixXd; 
using Eigen::VectorXd; 
using namespace Eigen;

VectorXd C; 
double  radius;

int iter = 0, back_iter = 0, mon_iter = 0;
int exit_condition1, exit_condition2;
int back_condition1, back_condition2, back_condition3;
int monitoring_condition1, monitoring_condition2, monitoring_condition3;

double FX, FV, FZ, FX_prev, FZ_prev;
double tk, qk, ck, tk_plus, qk_plus, ck_plus;
double alpha_Y, alpha_X;
double eta, del, rho_Y, rho_X;
int max_iter;

VectorXd X, Y, Z, V, X_prev, Y_prev;
VectorXd grad_X, grad_X_prev, grad_Y, grad_Y_prev;
VectorXd sk, rk;

double Objective_Function(const Vector2d& X){
    return 15.0*pow(X(0)*X(0)-1 , 2) + 1.0*pow(X(1)*X(1)-2 , 2) + 4*X(0)*X(1) + X(0) + X(1);
}

Vector2d Gradient(const Vector2d& X){
    Vector2d grad_X;
    grad_X << 4*X(1) + 60*X(0)*(X(0)*X(0) - 1) + 1,
              4*X(0) + 4*X(1)*(X(1)*X(1) - 2) + 1;
    
    return grad_X;
}


Vector2d Projection(const Vector2d& X){
    if((X-C).norm() > radius){
        return C + radius * (X-C)/(X-C).norm();
    }
    else{
        return X;
    }
}


void Projected_GD( Ref<VectorXd> X ){

    Y = X;
    
    FX = Objective_Function(X); FX_prev = FX;
    grad_X = Gradient(X);
    grad_Y = Gradient(Y);

    X_prev = X;
    Y_prev = VectorXd::Zero(Y.size());
    grad_Y_prev = VectorXd::Zero(Y.size());

    ////////////////////////////////////////////////////

    tk = 1;
    tk_plus = 1;
    qk = 1;
    ck = FX;
    eta = 0.4;
    del = 0.001;
    
    //iteration start
    max_iter = 100;
    iter = 0;
    while(1){
        iter = iter + 1;

        // backtracking
        grad_Y = Gradient(Y);
        sk = Y - Y_prev;
        rk = grad_Y - grad_Y_prev;
        alpha_Y = abs( (sk.transpose()*rk/(rk.transpose()*rk)).value() );

        rho_Y = 0.5;
        back_iter = 0;
        while (1){
            back_iter = back_iter + 1;
            Z = Projection(Y - alpha_Y * grad_Y);
            alpha_Y = rho_Y * alpha_Y;
            FZ = Objective_Function(Z);

            // Backtracking exit conditions
            back_condition1 = (ck - FZ) >= del * (Y-Z).transpose()*(Y-Z);
            back_condition2 = grad_Y.norm() < 1e-15;
            back_condition3 = back_iter > 10;
            if (back_condition1 || back_condition2 || back_condition3){
				break;
			} 
        }

        if (back_condition1){
            X = Z;
            FX = FZ;
        }
        else{
            grad_X = Gradient(X);
            sk = X - Y_prev;
            rk = grad_X - grad_Y_prev;
            alpha_X = abs( (sk.transpose()*rk/(rk.transpose()*rk)).value() );

            rho_X = 0.5;
            mon_iter = 0;
            while (1){
                mon_iter = mon_iter + 1;
                V = Projection(X - alpha_X * grad_X);
                alpha_X = rho_X * alpha_X;
                FV = Objective_Function(V);

                monitoring_condition1 = (ck - FV) >= del * (X-V).transpose()*(X-V);
                monitoring_condition2 = grad_X.norm() < 1e-15;
                monitoring_condition3 = mon_iter > 10;
                if (monitoring_condition1 || monitoring_condition2 || monitoring_condition3){
                    break;
                } 
            }

            if (FZ <= FV){
                X = Z;
                FX = FZ;
            }            
            else{
                X = V;
                FX = FV;
            }

        }

        tk_plus = (1+sqrt(1+4*tk*tk)) / 2.0;
        qk_plus = eta*qk + 1;
        ck_plus = (eta*qk*ck + FX)/qk_plus;

        exit_condition1 = pow(ck_plus-ck, 2) < 1e-6;
        exit_condition2 = iter >= max_iter; 
        if (exit_condition1 || exit_condition2){
            break;
        }

        // Nesterov step update
        Y_prev = Y;
        grad_Y_prev = grad_Y; 
        Y = X + tk/tk_plus*(Z-X) + (tk-1)/tk_plus*(X-X_prev);

        // Next step updates
        X_prev = X;
        tk = tk_plus;
        qk = qk_plus;
        ck = ck_plus;
    }

}

int main() {
    C = Vector2d(0, 1.2);
    radius = 0.5;
    
    X = Vector2d(-0.07, -1.45); 

    auto startTime = std::chrono::system_clock::now();
    Projected_GD(X);
    auto endTime = std::chrono::system_clock::now();
    auto time_measured = std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - startTime);

    std::cout << "time_measured: " << 1e-9 * time_measured.count() << std::endl;
    std::cout << "iter: " << iter << std::endl;
    std::cout << "X: " << X.transpose() << std::endl;
    std::cout << "objective: " << FX << std::endl;
    
    return 0;
}