function Fig10()

% Load data
M1 = cell2mat(struct2cell(load('./Data/M1.mat')));
M1_1 = M1(:,1:2501,1); M1_2 = M1(:,1:2501,2); 
M1_3 = M1(:,1:2501,3); M1_4 = M1(:,1:2501,4);
M1_5 = M1(:,1:2501,5); M1_6 = M1(:,1:2501,6); 
M1_7 = M1(:,1:2501,7); M1_8 = M1(:,1:2501,8);
time = 0:0.01:25;
M1_1_integrals = trapz(time, M1_1, 2); 
M1_2_integrals = trapz(time, M1_2, 2); 
M1_3_integrals = trapz(time, M1_3, 2); 
M1_4_integrals = trapz(time, M1_4, 2); 
M1_5_integrals = trapz(time, M1_5, 2); 
M1_6_integrals = trapz(time, M1_6, 2);
M1_7_integrals = trapz(time, M1_7, 2); 
M1_8_integrals = trapz(time, M1_8, 2); 
scale1 = 1e-7;
M1_mean = [mean(M1_1_integrals); mean(M1_2_integrals); mean(M1_3_integrals); ...
    mean(M1_4_integrals); mean(M1_5_integrals); mean(M1_6_integrals); ...
    mean(M1_7_integrals); mean(M1_8_integrals)].*scale1;

M3 = cell2mat(struct2cell(load('./Data/M3.mat')));
M3_1 = M3(:,1:2501,1); M3_2 = M3(:,1:2501,2); 
M3_3 = M3(:,1:2501,3); M3_4 = M3(:,1:2501,4);
M3_5 = M3(:,1:2501,5); M3_6 = M3(:,1:2501,6); 
M3_7 = M3(:,1:2501,7); M3_8 = M3(:,1:2501,8);
[~, col_M3_1] = max(M3_1, [], 2); [~, col_M3_2] = max(M3_2, [], 2);
[~, col_M3_3] = max(M3_3, [], 2); [~, col_M3_4] = max(M3_4, [], 2);
[~, col_M3_5] = max(M3_5, [], 2); [~, col_M3_6] = max(M3_6, [], 2);
[~, col_M3_7] = max(M3_7, [], 2); [~, col_M3_8] = max(M3_8, [], 2);
scale2 = 0.01;
M3_mean = [mean(col_M3_1); mean(col_M3_2); mean(col_M3_3); mean(col_M3_4); ...
    mean(col_M3_5); mean(col_M3_6); mean(col_M3_7); mean(col_M3_8)].*scale2;

Death = cell2mat(struct2cell(load('./Data/death_times.mat')));

% Figure
Color = ColorMatrix();
FontSize = 14;
LineWidth = 2;
MarkerSize = 9;

figure()
hold on; grid on; box on;
Colors = {cell2mat(Color(3,9)); cell2mat(Color(3,10)); cell2mat(Color(3,1)); ...
    cell2mat(Color(5,1)); cell2mat(Color(7,1)); cell2mat(Color(3,11)); ...
    cell2mat(Color(3,13)); cell2mat(Color(5,13))};
median_times = zeros(8, 1);
for i = 1:8
    idx = (i-1)*100 + (1:100);
    times = Death(idx);
    times = times(:);
    t_sorted = sort(times);
    n = length(t_sorted);
    survival = (n:-1:1)' / n; 
    x_data = [0; t_sorted];
    y_data = [1; survival];
    stairs(x_data, y_data*100, 'Color', Colors{i}, 'LineWidth', LineWidth);
    [~, idx_closest] = min(abs(y_data - 0.5));
    median_times(i) = x_data(idx_closest);
end
for i = 1:8
    scatter(median_times(i), 50, MarkerSize+120, ...
        'MarkerFaceColor', Colors{i}, 'MarkerEdgeColor', [1 1 1]);
end
legend('Mouse 1', 'Mouse 2', 'Mouse 3', 'Mouse 4', 'Mouse 5', 'Mouse 6', 'Mouse 7', 'Mouse 8', ...
    'Location', 'southwest', 'NumColumns', 1, 'FontSize', FontSize);
xlim([20 30]);
xlabel('Time (days)'); ylabel('Survival rate (%)');
set(gca, 'FontSize', FontSize, 'FontName', 'Arial');
set(gcf,'unit','centimeters','position',[3 10 20 20])
print('Figure/Fig10A','-dpng','-r200');

figure()
subplot(2, 2, 1);
hold on; grid on; box on;
for i = 1:8
    scatter(M1_mean(i), median_times(i), MarkerSize+100, 'MarkerEdgeColor', [1 1 1], ...
        'MarkerFaceColor', Colors{i}, 'LineWidth', 1.5);
end
p = polyfit(M1_mean, median_times, 1);
y_mean = polyval(p, M1_mean);
plot(M1_mean, y_mean, '-', 'Color', [0 0 0], 'LineWidth', 1.5);
[R_matrix, P_matrix] = corrcoef(M1_mean, median_times);
R_value = R_matrix(1,2);
P_value = P_matrix(1,2);
if P_value < 0.001
    p_text = 'p < 0.001';
elseif P_value < 0.01
    p_text = sprintf('p = %.2f', P_value);
else
    p_text = sprintf('p = %.2f', P_value);
end
text(0.1, 0.8, sprintf('R = %.2f\n%s', R_value, p_text), ...
     'Color', [0 0 0], 'FontSize', FontSize+2, 'Units', 'normalized');
legend('Mouse 1', 'Mouse 2', 'Mouse 3', 'Mouse 4', 'Mouse 5', 'Mouse 6', 'Mouse 7', 'Mouse 8', ...
    'Location', 'southeast', 'NumColumns', 1, 'FontSize', FontSize-3);
title('Group');
xlim([2.5 5.5]);
xlabel('M1-type macrophages AUC (\times 10^7 cell)'); ylabel('M-OS (days)');
set(gca, 'FontSize', FontSize, 'FontName', 'Arial');

subplot(2, 2, 2);
hold on; grid on; box on;
y_data = Death;
x_data = [M1_1_integrals; M1_2_integrals; M1_3_integrals; M1_4_integrals; ...
    M1_5_integrals; M1_6_integrals; M1_7_integrals; M1_8_integrals].*scale1;
hold on; grid on; box on;
M1_data = {M1_1_integrals, M1_2_integrals, M1_3_integrals, M1_4_integrals, ...
    M1_5_integrals, M1_6_integrals, M1_7_integrals, M1_8_integrals};
T_data = {Death(1:100), Death(101:200), Death(201:300), Death(301:400), ...
          Death(401:500), Death(501:600), Death(601:700), Death(701:800)};
for i = 1:length(M1_data)
    x = M1_data{i} .* scale1;
    y = T_data{i};
    scatter(x, y, MarkerSize+70, 'MarkerEdgeColor', [1 1 1], 'MarkerFaceColor', Colors{i});
end
p = polyfit(x_data, y_data, 1);
y_fit = polyval(p, x_data);
plot(x_data, y_fit, '-', 'Color', [0 0 0], 'LineWidth', LineWidth);
[R_matrix, P_matrix] = corrcoef(x_data, y_data);
R_value = R_matrix(1,2);
P_value = P_matrix(1,2);
if P_value < 0.001
    p_text = 'p < 0.001';
elseif P_value < 0.01
    p_text = sprintf('p = %.3f', P_value);
else
    p_text = sprintf('p = %.2f', P_value);
end
text(0.7, 0.2, sprintf('R = %.2f\n%s', R_value, p_text), ...
     'Color', [0 0 0], 'FontSize', FontSize+2, 'Units', 'normalized');
xlim([1.3 11]); ylim([21 30]);
title('Individual');
xlabel('M1-type macrophages AUC (\times 10^7 cell)'); ylabel('Death time (days)');
set(gca, 'FontSize', FontSize, 'FontName', 'Arial');

subplot(2, 2, 3);
hold on; grid on; box on;
for i = 1:8
    scatter(M3_mean(i), median_times(i), MarkerSize+100, 'MarkerEdgeColor', [1 1 1], ...
        'MarkerFaceColor', Colors{i}, 'LineWidth', 1.5);
end
p = polyfit(M3_mean, median_times, 1);
y_mean = polyval(p, M3_mean);
plot(M3_mean, y_mean, '-', 'Color', [0 0 0], 'LineWidth', 1.5);
[R_matrix, P_matrix] = corrcoef(M3_mean, median_times);
R_value = R_matrix(1,2);
P_value = P_matrix(1,2);
if P_value < 0.001
    p_text = 'p < 0.001';
elseif P_value < 0.01
    p_text = sprintf('p = %.2f', P_value);
else
    p_text = sprintf('p = %.2f', P_value);
end
text(0.1, 0.8, sprintf('R = %.2f\n%s', R_value, p_text), ...
     'Color', [0 0 0], 'FontSize', FontSize+2, 'Units', 'normalized');
legend('Mouse 1', 'Mouse 2', 'Mouse 3', 'Mouse 4', 'Mouse 5', 'Mouse 6', 'Mouse 7', 'Mouse 8', ...
    'Location', 'southeast', 'NumColumns', 1, 'FontSize', FontSize-3);
title('Group');
xlim([7 12]);
xlabel('Peak time of M3-type macrophages'); ylabel('M-OS (days)');
set(gca, 'FontSize', FontSize, 'FontName', 'Arial');

subplot(2, 2, 4);
hold on; grid on; box on;
y_data = Death;
x_data = [col_M3_1; col_M3_2; col_M3_3; col_M3_4; ...
    col_M3_5; col_M3_6; col_M3_7; col_M3_8].*scale2;
hold on; grid on; box on;
col_data = {col_M3_1, col_M3_2, col_M3_3, col_M3_4, col_M3_5, col_M3_6, col_M3_7, col_M3_8};
T_data = {Death(1:100), Death(101:200), Death(201:300), Death(301:400), ...
          Death(401:500), Death(501:600), Death(601:700), Death(701:800)};
for i = 1:length(col_data)
    x = col_data{i} .* scale2;
    y = T_data{i};
    scatter(x, y, MarkerSize+70, 'MarkerEdgeColor', [1 1 1], 'MarkerFaceColor', Colors{i});
end
p = polyfit(x_data, y_data, 1);
y_fit = polyval(p, x_data);
plot(x_data, y_fit, '-', 'Color', [0 0 0], 'LineWidth', 1.5);
[R_matrix, P_matrix] = corrcoef(x_data, y_data);
R_value = R_matrix(1,2);
P_value = P_matrix(1,2);
if P_value < 0.001
    p_text = 'p < 0.001';
elseif P_value < 0.01
    p_text = sprintf('p = %.3f', P_value);
else
    p_text = sprintf('p = %.2f', P_value);
end
text(0.7, 0.2, sprintf('R = %.2f\n%s', R_value, p_text), ...
     'Color', [0 0 0], 'FontSize', FontSize+2, 'Units', 'normalized');
xlabel('Peak time of M3-type macrophages');
ylabel('Death time (days)');
title('Individual');
ylim([21 30]); xlim([5 25]);
set(gca, 'FontSize', FontSize, 'FontName', 'Arial');


set(gcf,'unit','centimeters','position',[22 10 26 20])
print('Figure/Fig10BCDE','-dpng','-r200');
