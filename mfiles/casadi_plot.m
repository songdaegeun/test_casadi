function main
    global XX YY ZZ x_opt y_opt circle_x circle_y line_x line_y;
    load_all('rosenbrock_contour.mat', 'rosenbrock_optimal_solution.mat', ...
            'rosenbrock_constraint_1.mat', 'rosenbrock_constraint_2.mat')
    
    [x_opt, y_opt] = rosenbrock_result();
    
    disp(['Optimal x: ', num2str(x_opt)]);
    disp(['Optimal y: ', num2str(y_opt)]);
end