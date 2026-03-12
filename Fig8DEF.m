function Fig8DEF()

% close all; clc; clear

% Load data
Q = 2501;
T = cell2mat(struct2cell(load('./Data/Tumor.mat')));
T_1 = T(:,1:Q,1); T_2 = T(:,1:Q,2); T_3 = T(:,1:Q,3); T_4 = T(:,1:Q,4);
T_5 = T(:,1:Q,5); T_6 = T(:,1:Q,6); T_7 = T(:,1:Q,7); T_8 = T(:,1:Q,8);
time = 0:0.01:25;
T_1_integrals = trapz(time, T_1, 2); std_value_T_1 = std(T_1_integrals);
T_2_integrals = trapz(time, T_2, 2); std_value_T_2 = std(T_2_integrals);
T_3_integrals = trapz(time, T_3, 2); std_value_T_3 = std(T_3_integrals);
T_4_integrals = trapz(time, T_4, 2); std_value_T_4 = std(T_4_integrals);
T_5_integrals = trapz(time, T_5, 2); std_value_T_5 = std(T_5_integrals);
T_6_integrals = trapz(time, T_6, 2); std_value_T_6 = std(T_6_integrals);
T_7_integrals = trapz(time, T_7, 2); std_value_T_7 = std(T_7_integrals);
T_8_integrals = trapz(time, T_8, 2); std_value_T_8 = std(T_8_integrals);

M3 = cell2mat(struct2cell(load('./Data/M3.mat')));
M3_1 = M3(:,1:Q,1); M3_2 = M3(:,1:Q,2); M3_3 = M3(:,1:Q,3); M3_4 = M3(:,1:Q,4);
M3_5 = M3(:,1:Q,5); M3_6 = M3(:,1:Q,6); M3_7 = M3(:,1:Q,7); M3_8 = M3(:,1:Q,8);

[~, col_M3_1] = max(M3_1, [], 2); [~, col_M3_2] = max(M3_2, [], 2);
[~, col_M3_3] = max(M3_3, [], 2); [~, col_M3_4] = max(M3_4, [], 2);
[~, col_M3_5] = max(M3_5, [], 2); [~, col_M3_6] = max(M3_6, [], 2);
[~, col_M3_7] = max(M3_7, [], 2); [~, col_M3_8] = max(M3_8, [], 2);

scale1 = 0.01; scale2 = 1e-9;
[f_M3_1, xi_M3_1] = ksdensity(col_M3_1.*scale1);
[f_M3_2, xi_M3_2] = ksdensity(col_M3_2.*scale1);
[f_M3_3, xi_M3_3] = ksdensity(col_M3_3.*scale1);
[f_M3_4, xi_M3_4] = ksdensity(col_M3_4.*scale1);
[f_M3_5, xi_M3_5] = ksdensity(col_M3_5.*scale1);
[f_M3_6, xi_M3_6] = ksdensity(col_M3_6.*scale1);
[f_M3_7, xi_M3_7] = ksdensity(col_M3_7.*scale1);
[f_M3_8, xi_M3_8] = ksdensity(col_M3_8.*scale1);

figure()
Color = ColorMatrix();
LineWidth = 2;
FontSize = 14;
MarkerSize = 40;
Colors = {cell2mat(Color(3,9)); cell2mat(Color(3,10)); cell2mat(Color(3,1)); ...
    cell2mat(Color(5,1)); cell2mat(Color(7,1)); cell2mat(Color(3,11)); ...
    cell2mat(Color(3,13)); cell2mat(Color(5,13))};

subplot(1,3,1)
hold on; grid on; box on;
xi_data = {xi_M3_1, xi_M3_2, xi_M3_3, xi_M3_4, xi_M3_5, xi_M3_6, xi_M3_7, xi_M3_8};
f_data = {f_M3_1, f_M3_2, f_M3_3, f_M3_4, f_M3_5, f_M3_6, f_M3_7, f_M3_8};
for i = 1:length(xi_data)
    plot(xi_data{i}, f_data{i}, 'LineWidth', LineWidth, 'Color', Colors{i});
end
legend('Mouse 1', 'Mouse 2', 'Mouse 3', 'Mouse 4', 'Mouse 5', 'Mouse 6', 'Mouse 7', 'Mouse 8', ...
    'Location', 'northeast', 'NumColumns', 2, 'FontSize', FontSize-2);
xlabel('Peak time of M3-type macrophages');
ylabel('Probability density');
title('Distribution curve');
set(gca, 'FontSize', FontSize, 'FontName', 'Arial');

subplot(1,3,2)
T_mean = [mean(T_1_integrals); mean(T_2_integrals); mean(T_3_integrals); ...
    mean(T_4_integrals); mean(T_5_integrals); mean(T_6_integrals); ...
    mean(T_7_integrals); mean(T_8_integrals)].*scale2;
M3_mean = [mean(col_M3_1); mean(col_M3_2); mean(col_M3_3); mean(col_M3_4); ...
    mean(col_M3_5); mean(col_M3_6); mean(col_M3_7); mean(col_M3_8)].*scale1;
hold on; grid on; box on;
for i = 1:8
    scatter(M3_mean(i), T_mean(i), MarkerSize+100, 'MarkerEdgeColor', [1 1 1], ...
        'MarkerFaceColor', Colors{i}, 'LineWidth', 1.5);
end
p = polyfit(M3_mean, T_mean, 1);
y_mean = polyval(p, M3_mean);
plot(M3_mean, y_mean, '-', 'Color', [0 0 0], 'LineWidth', 1.5);
[R_matrix, P_matrix] = corrcoef(M3_mean, T_mean);
R_value = R_matrix(1,2);
P_value = P_matrix(1,2);
if P_value < 0.001
    p_text = 'p < 0.001';
elseif P_value < 0.01
    p_text = sprintf('p = %.3f', P_value);
else
    p_text = sprintf('p = %.2f', P_value);
end
text(0.1, 0.2, sprintf('R = %.2f\n%s', R_value, p_text), ...
     'Color', [0 0 0], 'FontSize', FontSize+2, 'Units', 'normalized');
xlim([7 11.5]); ylim([1.25 3]);
xlabel('Peak time of M3-type macrophages');
ylabel('Tumor AUC (\times 10^9 mm^3)');
legend('Mouse 1', 'Mouse 2', 'Mouse 3', 'Mouse 4', 'Mouse 5', 'Mouse 6', 'Mouse 7', 'Mouse 8', ...
    'Location', 'northeast', 'NumColumns', 2, 'FontSize', FontSize-2);
set(gca, 'FontSize', FontSize, 'FontName', 'Arial');
title('Group');

subplot(1,3,3)
y_data = [T_1_integrals; T_2_integrals; T_3_integrals; T_4_integrals; ...
    T_5_integrals; T_6_integrals; T_7_integrals; T_8_integrals].*scale2;
x_data = [col_M3_1; col_M3_2; col_M3_3; col_M3_4; ...
    col_M3_5; col_M3_6; col_M3_7; col_M3_8].*scale1;
hold on; grid on; box on;
col_data = {col_M3_1, col_M3_2, col_M3_3, col_M3_4, col_M3_5, col_M3_6, col_M3_7, col_M3_8};
T_data = {T_1_integrals, T_2_integrals, T_3_integrals, T_4_integrals, ...
          T_5_integrals, T_6_integrals, T_7_integrals, T_8_integrals};
for i = 1:length(col_data)
    x = col_data{i} .* scale1;
    y = T_data{i} .* scale2;
    scatter(x, y, MarkerSize, 'MarkerEdgeColor', [1 1 1], 'MarkerFaceColor', Colors{i});
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
text(0.7, 0.8, sprintf('R = %.2f\n%s', R_value, p_text), ...
     'Color', [0 0 0], 'FontSize', FontSize+2, 'Units', 'normalized');
xlabel('Peak time of M3-type macrophages');
ylabel('Tumor AUC (\times 10^9 mm^3)');
title('Individual');
ylim([1.2 3.5]);
set(gca, 'FontSize', FontSize, 'FontName', 'Arial');

set(gcf,'unit','centimeters','position',[3 10 45 10])
print('Figure/Fig8DEF','-dpng','-r200')
