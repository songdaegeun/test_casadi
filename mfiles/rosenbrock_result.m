function [x_opt, y_opt] = rosenbrock_result()
    global XX YY ZZ x_opt y_opt circle_x circle_y line_x line_y;

    % Rosenbrock contour plot
    % Load precomputed data from CasADi C++ output
    figure;
    contour(XX, YY, ZZ, 100);
    colormap('parula');
    colorbar;
    hold on;
    plot(x_opt, y_opt, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    plot(circle_x, circle_y, 'r', 'LineWidth', 2);
    plot(line_x, line_y, 'r', 'LineWidth', 2);
    
    xlabel('x'); ylabel('y'); title('Rosenbrock Contour');
    xlim([0, 1.5]);
    ylim([-0.5, 1.5]);
    grid on;
end
