function main
    [x_opt, y_opt] = rosenbrock_result('rosenbrock_data.mat')
    disp(['Optimal x: ', num2str(x_opt)]);
    disp(['Optimal y: ', num2str(y_opt)]);
end