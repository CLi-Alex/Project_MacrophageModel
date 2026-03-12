function Fig9EFGHIJKL()

% Load data
Matrix = cell2mat(struct2cell(load('./Data/Matrix.mat')));
Alpha = [Matrix(1,:,1) Matrix(1,:,2) Matrix(1,:,3) Matrix(1,:,4) ...
    Matrix(1,:,5) Matrix(1,:,6) Matrix(1,:,7) Matrix(1,:,8)]';
K = [Matrix(2,:,1) Matrix(2,:,2) Matrix(2,:,3) Matrix(2,:,4) ...
    Matrix(2,:,5) Matrix(2,:,6) Matrix(2,:,7) Matrix(2,:,8)]';
M0 = [Matrix(3,:,1) Matrix(3,:,2) Matrix(3,:,3) Matrix(3,:,4) ...
    Matrix(3,:,5) Matrix(3,:,6) Matrix(3,:,7) Matrix(3,:,8)]';
Kappa = [Matrix(4,:,1) Matrix(4,:,2) Matrix(4,:,3) Matrix(4,:,4) ...
    Matrix(4,:,5) Matrix(4,:,6) Matrix(4,:,7) Matrix(4,:,8)]';
Beta_13 = [Matrix(5,:,1) Matrix(5,:,2) Matrix(5,:,3) Matrix(5,:,4) ...
    Matrix(5,:,5) Matrix(5,:,6) Matrix(5,:,7) Matrix(5,:,8)]';
Beta_32 = [Matrix(6,:,1) Matrix(6,:,2) Matrix(6,:,3) Matrix(6,:,4) ...
    Matrix(6,:,5) Matrix(6,:,6) Matrix(6,:,7) Matrix(6,:,8)]';
Eta = [Matrix(7,:,1) Matrix(7,:,2) Matrix(7,:,3) Matrix(7,:,4) ...
    Matrix(7,:,5) Matrix(7,:,6) Matrix(7,:,7) Matrix(7,:,8)]';
P = [Matrix(8,:,1) Matrix(8,:,2) Matrix(8,:,3) Matrix(8,:,4) ...
    Matrix(8,:,5) Matrix(8,:,6) Matrix(8,:,7) Matrix(8,:,8)]';

Death = cell2mat(struct2cell(load('./Data/death_times.mat')));

% Figure
figure()
Color = ColorMatrix();
FontSize = 14;
LineWidth = 2;
MarkerSize = 9;

subplot(2, 4, 1);
hold on; grid on; box on;
Q = median(Alpha);
low_idx = Alpha < Q;
high_idx = ~low_idx;
[~, t_low, s_low] = ecdf(Death(low_idx), 'Function', 'survivor');
[~, t_high, s_high] = ecdf(Death(high_idx), 'Function', 'survivor');
[~, idx_low] = min(abs(s_low - 0.5));
[~, idx_high] = min(abs(s_high - 0.5));
stairs(t_low, s_low*100, 'b-', 'LineWidth', 2);
stairs(t_high, s_high*100, 'r-', 'LineWidth', 2);
scatter(t_low(idx_low), 50, MarkerSize+100, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', [1 1 1]);
scatter(t_high(idx_high), 50, MarkerSize+100, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', [1 1 1]);
text(25.3, 57, sprintf('%.1f days', t_high(idx_low)), 'FontSize', FontSize, 'FontWeight', 'bold', 'Color', 'b');
text(20.2, 44, sprintf('%.1f days', t_high(idx_high)), 'FontSize', FontSize, 'FontWeight', 'bold', 'Color', 'r');
title('\alpha'); xlabel('Time (days)'); ylabel('Survival rate（%）');
xlim([20 30]);
legend('Low', 'High', 'Location', 'northeast');
set(gca, 'FontSize', FontSize, 'FontName', 'Arial');

subplot(2, 4, 2);
hold on; grid on; box on;
Q = median(K);
low_idx = K < Q;
high_idx = ~low_idx;
[~, t_low, s_low] = ecdf(Death(low_idx), 'Function', 'survivor');
[~, t_high, s_high] = ecdf(Death(high_idx), 'Function', 'survivor');
[~, idx_low] = min(abs(s_low - 0.5));
[~, idx_high] = min(abs(s_high - 0.5));
stairs(t_low, s_low*100, 'b-', 'LineWidth', 2);
stairs(t_high, s_high*100, 'r-', 'LineWidth', 2);
scatter(t_low(idx_low), 50, MarkerSize+100, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', [1 1 1]);
scatter(t_high(idx_high), 50, MarkerSize+100, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', [1 1 1]);
text(20.5, 44, sprintf('%.1f days', t_high(idx_low)), 'FontSize', FontSize, 'FontWeight', 'bold', 'Color', 'b');
text(25.1, 57, sprintf('%.1f days', t_high(idx_high)), 'FontSize', FontSize, 'FontWeight', 'bold', 'Color', 'r');
title('K'); xlabel('Time (days)'); ylabel('Survival rate（%）');
xlim([20 30]);
legend('Low', 'High', 'Location', 'northeast');
set(gca, 'FontSize', FontSize, 'FontName', 'Arial');

subplot(2, 4, 3);
hold on; grid on; box on;
Q = median(M0);
low_idx = M0 < Q;
high_idx = ~low_idx;
[~, t_low, s_low] = ecdf(Death(low_idx), 'Function', 'survivor');
[~, t_high, s_high] = ecdf(Death(high_idx), 'Function', 'survivor');
[~, idx_low] = min(abs(s_low - 0.5));
[~, idx_high] = min(abs(s_high - 0.5));
stairs(t_low, s_low*100, 'b-', 'LineWidth', 2);
stairs(t_high, s_high*100, 'r-', 'LineWidth', 2);
scatter(t_low(idx_low), 50, MarkerSize+100, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', [1 1 1]);
scatter(t_high(idx_high), 50, MarkerSize+100, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', [1 1 1]);
text(20.5, 44, sprintf('%.1f days', t_high(idx_low)), 'FontSize', FontSize, 'FontWeight', 'bold', 'Color', 'b');
text(25.1, 57, sprintf('%.1f days', t_high(idx_high)), 'FontSize', FontSize, 'FontWeight', 'bold', 'Color', 'r');
title('M_0'); xlabel('Time (days)'); ylabel('Survival rate（%）');
xlim([20 30]);
legend('Low', 'High', 'Location', 'northeast');
set(gca, 'FontSize', FontSize, 'FontName', 'Arial');

subplot(2, 4, 4);
hold on; grid on; box on;
Q = median(Kappa);
low_idx = Kappa < Q;
high_idx = ~low_idx;
[~, t_low, s_low] = ecdf(Death(low_idx), 'Function', 'survivor');
[~, t_high, s_high] = ecdf(Death(high_idx), 'Function', 'survivor');
[~, idx_low] = min(abs(s_low - 0.5));
[~, idx_high] = min(abs(s_high - 0.5));
stairs(t_low, s_low*100, 'b-', 'LineWidth', 2);
stairs(t_high, s_high*100, 'r-', 'LineWidth', 2);
scatter(t_low(idx_low), 50, MarkerSize+100, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', [1 1 1]);
scatter(t_high(idx_high), 50, MarkerSize+100, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', [1 1 1]);
text(20.5, 44, sprintf('%.1f days', t_high(idx_low)), 'FontSize', FontSize, 'FontWeight', 'bold', 'Color', 'b');
text(25.1, 57, sprintf('%.1f days', t_high(idx_high)), 'FontSize', FontSize, 'FontWeight', 'bold', 'Color', 'r');
title('\kappa'); xlabel('Time (days)'); ylabel('Survival rate（%）');
xlim([20 30]);
legend('Low', 'High', 'Location', 'northeast');
set(gca, 'FontSize', FontSize, 'FontName', 'Arial');

subplot(2, 4, 5);
hold on; grid on; box on;
Q = median(Beta_13);
low_idx = Beta_13 < Q;
high_idx = ~low_idx;
[~, t_low, s_low] = ecdf(Death(low_idx), 'Function', 'survivor');
[~, t_high, s_high] = ecdf(Death(high_idx), 'Function', 'survivor');
[~, idx_low] = min(abs(s_low - 0.5));
[~, idx_high] = min(abs(s_high - 0.5));
stairs(t_low, s_low*100, 'b-', 'LineWidth', 2);
stairs(t_high, s_high*100, 'r-', 'LineWidth', 2);
scatter(t_low(idx_low), 50, MarkerSize+100, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', [1 1 1]);
scatter(t_high(idx_high), 50, MarkerSize+100, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', [1 1 1]);
text(25.1, 57, sprintf('%.1f days', t_high(idx_low)), 'FontSize', FontSize, 'FontWeight', 'bold', 'Color', 'b');
text(20.5, 44, sprintf('%.1f days', t_high(idx_high)), 'FontSize', FontSize, 'FontWeight', 'bold', 'Color', 'r');
title('\beta_{13}'); xlabel('Time (days)'); ylabel('Survival rate（%）');
xlim([20 30]);
legend('Low', 'High', 'Location', 'northeast');
set(gca, 'FontSize', FontSize, 'FontName', 'Arial');

subplot(2, 4, 6);
hold on; grid on; box on;
Q = median(Beta_32);
low_idx = Beta_32 < Q;
high_idx = ~low_idx;
[~, t_low, s_low] = ecdf(Death(low_idx), 'Function', 'survivor');
[~, t_high, s_high] = ecdf(Death(high_idx), 'Function', 'survivor');
[~, idx_low] = min(abs(s_low - 0.5));
[~, idx_high] = min(abs(s_high - 0.5));
stairs(t_low, s_low*100, 'b-', 'LineWidth', 2);
stairs(t_high, s_high*100, 'r-', 'LineWidth', 2);
scatter(t_low(idx_low), 50, MarkerSize+100, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', [1 1 1]);
scatter(t_high(idx_high), 50, MarkerSize+100, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', [1 1 1]);
text(25.1, 57, sprintf('%.1f days', t_high(idx_low)), 'FontSize', FontSize, 'FontWeight', 'bold', 'Color', 'b');
text(20.5, 44, sprintf('%.1f days', t_high(idx_high)), 'FontSize', FontSize, 'FontWeight', 'bold', 'Color', 'r');
title('\beta_{32}'); xlabel('Time (days)'); ylabel('Survival rate（%）');
xlim([20 30]);
legend('Low', 'High', 'Location', 'northeast');
set(gca, 'FontSize', FontSize, 'FontName', 'Arial');

subplot(2, 4, 7);
hold on; grid on; box on;
Q = median(Eta);
low_idx = Eta < Q;
high_idx = ~low_idx;
[~, t_low, s_low] = ecdf(Death(low_idx), 'Function', 'survivor');
[~, t_high, s_high] = ecdf(Death(high_idx), 'Function', 'survivor');
[~, idx_low] = min(abs(s_low - 0.5));
[~, idx_high] = min(abs(s_high - 0.5));
stairs(t_low, s_low*100, 'b-', 'LineWidth', 2);
stairs(t_high, s_high*100, 'r-', 'LineWidth', 2);
scatter(t_low(idx_low), 50, MarkerSize+100, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', [1 1 1]);
scatter(t_high(idx_high), 50, MarkerSize+100, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', [1 1 1]);
text(20.5, 44, sprintf('%.1f days', t_high(idx_low)), 'FontSize', FontSize, 'FontWeight', 'bold', 'Color', 'b');
text(25.1, 57, sprintf('%.1f days', t_high(idx_high)), 'FontSize', FontSize, 'FontWeight', 'bold', 'Color', 'r');
title('\eta'); xlabel('Time (days)'); ylabel('Survival rate（%）');
xlim([20 30]);
legend('Low', 'High', 'Location', 'northeast');
set(gca, 'FontSize', FontSize, 'FontName', 'Arial');

subplot(2, 4, 8);
hold on; grid on; box on;
Q = median(P);
low_idx = P < Q;
high_idx = ~low_idx;
[~, t_low, s_low] = ecdf(Death(low_idx), 'Function', 'survivor');
[~, t_high, s_high] = ecdf(Death(high_idx), 'Function', 'survivor');
[~, idx_low] = min(abs(s_low - 0.5));
[~, idx_high] = min(abs(s_high - 0.5));
stairs(t_low, s_low*100, 'b-', 'LineWidth', 2);
stairs(t_high, s_high*100, 'r-', 'LineWidth', 2);
scatter(t_low(idx_low), 50, MarkerSize+100, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', [1 1 1]);
scatter(t_high(idx_high), 50, MarkerSize+100, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', [1 1 1]);
text(25.1, 57, sprintf('%.1f days', t_high(idx_low)), 'FontSize', FontSize, 'FontWeight', 'bold', 'Color', 'b');
text(20.5, 44, sprintf('%.1f days', t_high(idx_high)), 'FontSize', FontSize, 'FontWeight', 'bold', 'Color', 'r');
title('p'); xlabel('Time (days)'); ylabel('Survival rate（%）');
xlim([20 30]);
legend('Low', 'High', 'Location', 'northeast');
set(gca, 'FontSize', FontSize, 'FontName', 'Arial');

set(gcf,'unit','centimeters','position',[3 10 50 18])
print('Figure/Fig9EFGHIJKL','-dpng','-r200');
