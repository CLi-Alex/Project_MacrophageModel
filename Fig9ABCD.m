function Fig9ABCD()

% Experimental data
Death = [100*ones(1,length(0:0.01:20)) ...
    87.5*ones(1,length(20.01:0.01:23)) ...
    62.5*ones(1,length(23.01:0.01:26)) ...
    12.5*ones(1,length(26.01:0.01:28)) ...
    0*ones(1,length(28.01:0.01:31))];

% Load data
Matrix = cell2mat(struct2cell(load('./Data/Matrix.mat')));
T = cell2mat(struct2cell(load('./Data/Tumor.mat')));
T_1 = zeros(100,3101);
T_2 = zeros(100,3101);
T_3 = zeros(100,3101);
T_4 = zeros(100,3101);
T_5 = zeros(100,3101);
T_6 = zeros(100,3101);
T_7 = zeros(100,3101);
T_8 = zeros(100,3101);
for j = 1:100
    T_1(j,:) = T(j,:,1)/Matrix(2,j,1);
    T_2(j,:) = T(j,:,2)/Matrix(2,j,2);
    T_3(j,:) = T(j,:,3)/Matrix(2,j,3);
    T_4(j,:) = T(j,:,4)/Matrix(2,j,4);
    T_5(j,:) = T(j,:,5)/Matrix(2,j,5);
    T_6(j,:) = T(j,:,6)/Matrix(2,j,6);
    T_7(j,:) = T(j,:,7)/Matrix(2,j,7);
    T_8(j,:) = T(j,:,8)/Matrix(2,j,8);
end


Time = 0:0.01:31;
rng(1); a = 0.025;
T_all = [T_1;T_2;T_3;T_4;T_5;T_6;T_7;T_8] + (a * randn(800,3101));

% PSO parameters
pso_params = struct('n_particles', 20, 'max_iter', 20, 'w', 0.7, 'c1', 2, 'c2', 2);
Q_range = [0.9 1];

% Run PSO optimization
[Q_opt, best_fitness, history] = pso_optimize_Q(T_all, Death, Time, Q_range, pso_params);
fprintf('Result：Q = %.4f, RMSE = %.4f\n', Q_opt, best_fitness);

death_times = calculate_death_times(T_all, Time, Q_opt);
simulated_survival = calculate_survival_curve(death_times, Time);

plot_results(Time, Death, simulated_survival, death_times, ...
    T_all, Q_opt, Q_range, history);

save('Data/death_times.mat', 'death_times');

end

function death_times = calculate_death_times(T_mean, Time, Q)
[n_samples, ~] = size(T_mean);
death_times = ones(n_samples, 1) * Time(end);

for i = 1:n_samples
    idx = find(T_mean(i, :) > Q, 1);
    if ~isempty(idx)
        death_times(i) = Time(idx);
    end
end
end

function survival_curve = calculate_survival_curve(death_times, Time)
n_samples = length(death_times);
survival_curve = zeros(1, length(Time));

for t_idx = 1:length(Time)
    alive_count = sum(death_times > Time(t_idx));
    survival_curve(t_idx) = (alive_count / n_samples) * 100;
end
end

function rmse = objective_function(Q, T_mean, Death, Time)
death_times = calculate_death_times(T_mean, Time, Q);
simulated_survival = calculate_survival_curve(death_times, Time);
rmse = sqrt(mean((simulated_survival - Death).^2));
end

function [Q_opt, best_fitness, history] = pso_optimize_Q(...
    T_mean, Death, Time, Q_range, params)

n_particles = params.n_particles;
max_iter = params.max_iter;
w = params.w;
c1 = params.c1;
c2 = params.c2;

particles.pos = Q_range(1) + rand(n_particles, 1) * (Q_range(2) - Q_range(1));
particles.vel = zeros(n_particles, 1);
particles.pbest = particles.pos;
particles.pbest_fit = inf(n_particles, 1);

gbest_fit = inf;
gbest_pos = particles.pos(1);
history = zeros(max_iter, 1);

for iter = 1:max_iter
    w_current = w * (1 - 0.5 * iter/max_iter);
    for i = 1:n_particles
        fitness = objective_function(particles.pos(i), T_mean, Death, Time);
        if fitness < particles.pbest_fit(i)
            particles.pbest_fit(i) = fitness;
            particles.pbest(i) = particles.pos(i);
        end
        if fitness < gbest_fit
            gbest_fit = fitness;
            gbest_pos = particles.pos(i);
        end
    end
    for i = 1:n_particles
        r1 = rand();
        r2 = rand();
        particles.vel(i) = w_current * particles.vel(i) + ...
            c1 * r1 * (particles.pbest(i) - particles.pos(i)) + ...
            c2 * r2 * (gbest_pos - particles.pos(i));
        particles.pos(i) = particles.pos(i) + particles.vel(i);
        if particles.pos(i) < Q_range(1)
            particles.pos(i) = Q_range(1);
            particles.vel(i) = -0.5 * particles.vel(i);
        elseif particles.pos(i) > Q_range(2)
            particles.pos(i) = Q_range(2);
            particles.vel(i) = -0.5 * particles.vel(i);
        end
    end
    history(iter) = gbest_fit;
    if mod(iter, 5) == 0
        fprintf('iteration %3d: RMSE = %.4f, Q = %.4f\n', iter, gbest_fit, gbest_pos);
    end
end

Q_opt = gbest_pos;
best_fitness = gbest_fit;
end

function plot_results(Time, Death, simulated_survival, death_times, T_mean, Q_opt, Q_range, history)

figure();
Color = ColorMatrix();
FontSize = 14;
LineWidth = 2;
MarkerSize = 9;

subplot(1, 4, 1);
hold on; grid on; box on;
plot(Time, Death, 'r-', 'LineWidth', LineWidth);
plot(Time, simulated_survival, 'b-.', 'LineWidth', LineWidth);
scatter(median(death_times), 50, MarkerSize+100, ...
    'MarkerFaceColor', 'b', 'MarkerEdgeColor', [1 1 1]);
text(15.75,50,sprintf('M-OS = %.1f days', median(death_times)), ...
    'FontSize', FontSize, 'FontWeight', 'bold', 'Color', 'b');
xlim([15 30]);
xlabel('Time (days)');
ylabel('Survival rate（%）');
title(sprintf('Survival curve (Q=%.4f)', Q_opt));
legend_labels = {
    ['Experimental data' newline '         (n = 8)']
    ['Model simulation' newline '       (n = 800)']
    };
legend(legend_labels, 'Location', 'southwest', 'FontSize', FontSize-3);
set(gca, 'FontSize', FontSize, 'FontName', 'Arial');

subplot(1, 4, 2);
hold on; grid on; box on;
plot(1:length(history), history, '-o', 'Color', cell2mat(Color(1,13)), ...
    'LineWidth', LineWidth, 'MarkerIndices', 1:1:length(history), ...
    'MarkerSize', MarkerSize, 'MarkerFaceColor', cell2mat(Color(1,13)), 'MarkerEdgeColor', [1 1 1]);
xlim([1 20]);
xlabel('Iterations');
ylabel('RMSE');
title('PSO convergence curve');
set(gca, 'FontSize', FontSize, 'FontName', 'Arial');

subplot(1, 4, 3);
hold on; grid on; box on;
plot(Time, T_mean(1:800, :), 'LineWidth', LineWidth-1.75);
yline(Q_opt, 'r:', 'LineWidth', LineWidth+1);
text(5, 0.85, sprintf('Q=%.4f', Q_opt), 'FontSize', FontSize+2, ...
    'FontWeight', 'bold', 'Color', 'r');
xlim([0 30]);
xlabel('Time (days)');
ylabel('Death risk');
title('Individual death risk');
set(gca, 'FontSize', FontSize, 'FontName', 'Arial');

subplot(1, 4, 4);
hold on; grid on; box on;
Q_test = linspace(Q_range(1), Q_range(2), 15);
rmse_test = zeros(size(Q_test));
for i = 1:length(Q_test)
    death_times_tmp = calculate_death_times(T_mean, Time, Q_test(i));
    survival_tmp = calculate_survival_curve(death_times_tmp, Time);
    rmse_test(i) = sqrt(mean((survival_tmp - Death).^2));
end
plot(Q_test, rmse_test, '-o', 'Color', 'b', 'LineWidth', LineWidth, ...
    'MarkerIndices', 1:1:length(Q_test), 'MarkerSize', MarkerSize, ...
    'MarkerFaceColor', 'b', 'MarkerEdgeColor', [1 1 1]);
xline(Q_opt, 'r:', 'LineWidth', LineWidth+1);
text(0.92, 7.25, sprintf('Q=%.4f', Q_opt), 'FontSize', FontSize+2, ...
    'FontWeight', 'bold', 'Color', 'r');
ylim([6 14]);
xlabel('Threshold Q');
ylabel('RMSE');
title('Q value sensitivity analysis');
set(gca, 'FontSize', FontSize, 'FontName', 'Arial');

set(gcf,'unit','centimeters','position',[3 10 50 9])
print('Figure/Fig9ABCD','-dpng','-r200');

end
