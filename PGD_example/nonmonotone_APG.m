%%
clear; clc;

%% Function Defining
x         = sym("x",[2 1]); 
% objective = 15*(x(1)^2-1)^2 + 1*(x(2)^2-2)^2 + 4*x(1)*x(2) + x(1) + x(2);
objective = 15*(x(1)^2-1)^2 + 1*(x(2)^2-2)^2 + 4*x(1)*x(2) + x(1) + x(2);
grad      = gradient(objective);

% constraint
global C radius;
C = [0; 1.2]; 
radius = 0.5;

%% Search
% 초기 시작 지점을 0으로 잡으면 안된다. rk가 0이 돼서 스텝 사이즈가 NaN이 됨.
% X = [0.1; 1.3];
% X = [0.025; 0.15];
% X = [0.1; 0.1];
% X = [-0.07; -1.45]; % eta=0.2,0.5,0.8
% X = [0.3; -0.6];
% X = [-0.01; -1.49];
X = [0.2; 1];
Y = X;

FX = double( subs(objective,x,X) ); F_data = FX;
FX_prev = FX;
grad_X = double( subs(grad,x,X) );
grad_Y = double( subs(grad,x,Y) );
X_prev = X;
Y_prev = [0; 0];
grad_Y_prev = [0; 0];
x_data(:,1) = X;

% tune eta<1, rho<1, del>0
tk = 1;
tk_plus = 1;
qk = 1;
eta = 0.4; % nonmonotonicity
ck = FX; c_data(1) = ck;
del = 0.001; %10;

tic;

% iterations
MAX_ITER = 100;
iter = 0;
while 1
    iter = iter +1

    % backtracking
    grad_Y = double( subs(grad,x,Y) ); % mag = norm(grad_y);
    sk = Y - Y_prev;
    rk = grad_Y - grad_Y_prev;
    % alpha_Y = abs(sk.'*sk/(sk.'*rk)); % long
    alpha_Y = abs(sk.'*rk/(rk.'*rk)); % short

    rho_Y = 0.5;
    back_iter = 0;
    while 1
        back_iter = back_iter + 1
        % Z = Y - alpha_Y * grad_Y;
        Z = projection(Y - alpha_Y * grad_Y);
        alpha_Y = rho_Y * alpha_Y;
        FZ = double( subs(objective,x,Z) );
        
        % Backtracking exit conditions
        back_condition1 = (ck-FZ >= del*(Y-Z).'*(Y-Z));
        back_condition2 = back_iter > 10;
        if (back_condition1 | back_condition2)
            break;
        end
    end

    if (back_condition1)
        X = Z;
        FX = FZ;

    else
        disp("MonitorMonitorMonitorMonitorMonitor")
        grad_X = double( subs(grad,x,X) );
        sk = X - Y_prev;
        rk = grad_X - grad_Y_prev;
        % alpha_X = abs(sk.'*sk/(sk.'*rk)); % long
        alpha_X = abs(sk.'*rk/(rk.'*rk)); % short

        rho_X = 0.5;
        mon_iter = 0;
        while 1
            mon_iter = mon_iter + 1
            % V = X - alpha_X * grad_X;
            V = projection(X - alpha_X * grad_X);
            alpha_X = rho_X * alpha_X;
            FV = double( subs(objective,x,V) );

            monitoring_condition1 = ck - FV >= del*(X-V).'*(X-V);
            monitoring_condition2 = mon_iter > 10;
            if (monitoring_condition1 | monitoring_condition2)
                break;
            end
        end

        if FZ <= FV
            X = Z;
            FX = FZ;
        else
            X = V;
            FX = FV;
        end

    end

    tk_plus = (1+sqrt(1+4*tk^2))/2;
    qk_plus = eta*qk + 1;
    ck_plus = (eta*qk*ck + FX)/qk_plus;

    x_data(:,iter+1) = X;
    F_data(iter+1) = FX;
    c_data(iter+1) = ck_plus;

    % exit. X,F가 정해졌으니 끝내는 부분이 바로 와줘야 함.
    % exit_condition1 = (FX-FX_prev)^2 < 1e-6;
    exit_condition1 = (ck_plus-ck)^2 < 1e-6;
    exit_condition2 = iter >= MAX_ITER;
    if (exit_condition1 | exit_condition2)
        break;
    end

    % Nesterov step update
    Y_prev = Y;
    grad_Y_prev = grad_Y;    
    Y = X + tk/tk_plus*(Z-X) + (tk-1)/tk_plus*(X-X_prev);

    % Next step updates
    X_prev = X;
    FX_prev = FX;
    tk = tk_plus;
    qk = qk_plus;
    ck = ck_plus;

end

toc;

%% plot
% figure; plot([0:iter],F_data); grid; xlabel("iterations"); ylabel("objective");
% figure; plot([0:iter],c_data); grid; xlabel("iterations"); ylabel("Avg. objective");
figure; hold on; grid; 
plot([0:iter],F_data, 'r',"LineWidth",1);
plot([0:iter],c_data, 'b',"LineWidth",1); 
% title("Objective vs. Combination");
legend("F(x_k)","C_k");
xlabel("iterations"); ylabel("objective Function");

f=figure; hold on; grid; xlabel("x_1"); ylabel("x_2"); %daspect([1 1 1]);
fcontour(objective,[-1.25 1.25, -2 2],"Levelstep",0.4); xlim([-1.25 1.25]);
plot(x_data(1,:),x_data(2,:), "k-o", "LineWidth",1);
% initial & last points
plot(x_data(1,1),x_data(2,1), "r-o", "LineWidth",1.5);
plot(x_data(1,end),x_data(2,end), "r-o", "LineWidth",1.5);
% constraint
fimplicit(@(x1,x2) (x1-C(1)).^2 + (x2-C(2)).^2 - radius^2, "r", "LineWidth",1);

% set(f,'Position',[100 100 600 200]);
% xlim([-0.25 1.22]); ylim([-2 -1]);
%% 해상도 설정하여 저장
% set(gcf, 'Renderer', 'Painters'); % 렌더러를 Painters로 설정
% exportgraphics(gcf,'asdf.png','Resolution',1200)

%% functions
function projected_X = projection(X)
    global C radius;

    if(norm(X-C) > radius)
        projected_X = C + radius * (X-C)/norm(X-C); 
    else
        projected_X = X;
    end

end
       