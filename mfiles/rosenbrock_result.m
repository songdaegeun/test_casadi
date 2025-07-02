function [x_opt, y_opt] = rosenbrock_result(dataFile)
    % Rosenbrock contour plot
    % Load precomputed data from CasADi C++ output
    load(dataFile);  % Loads XX, YY, ZZ, x_opt, y_opt
    figure;
    contour(XX, YY, ZZ, 100);
    colormap('parula');
    colorbar;
    hold on;
    plot(x_opt, y_opt, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    xlabel('x'); ylabel('y'); title('Rosenbrock Contour');
    grid on;
end
