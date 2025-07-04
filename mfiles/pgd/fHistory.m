% fHistory.m

% JSON 파일 읽기
jsonText = fileread('x0_x_opt.json');
data = jsondecode(jsonText);

num_horizon = 10;  

% figure 초기화
figure;
for h = 1:num_horizon
    subplot(2, 5, h);
    hold on;
    grid on;
    xlabel('Step');
    ylabel('Cost');
    title(sprintf('Horizon step %d', h-1));
end

% 데이터 누적 플롯
for t = 1:length(data)
    time_step = data(t).time_step;
    horizon_steps = data(t).horizon_steps;

    for h = 1:length(horizon_steps)
        step_data = horizon_steps(h);
        cost = step_data.cost;
        initial_x = step_data.initial_x;
        optimal_x = step_data.optimal_x;

        steps = 1:length(cost);

        % horizon 별 subplot에 추가 플롯
        subplot(2, 5, h);
        plot(steps, cost, '-', 'DisplayName', sprintf('t=%d', time_step));

        % initial_x(1), optimal_x(1)
        title(sprintf('Step %d\n', h-1));
    end
end
for h = 1:10
    subplot(2, 5, h);
    legend('show');
end