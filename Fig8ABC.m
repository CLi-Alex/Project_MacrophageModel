function Fig8ABC()

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

M1 = cell2mat(struct2cell(load('./Data/M1.mat')));
M1_1 = M1(:,1:Q,1); M1_2 = M1(:,1:Q,2); M1_3 = M1(:,1:Q,3); M1_4 = M1(:,1:Q,4);
M1_5 = M1(:,1:Q,5); M1_6 = M1(:,1:Q,6); M1_7 = M1(:,1:Q,7); M1_8 = M1(:,1:Q,8);
time = 0:0.01:25;
M1_1_integrals = trapz(time, M1_1, 2); 
M1_2_integrals = trapz(time, M1_2, 2); 
M1_3_integrals = trapz(time, M1_3, 2); 
M1_4_integrals = trapz(time, M1_4, 2); 
M1_5_integrals = trapz(time, M1_5, 2); 
M1_6_integrals = trapz(time, M1_6, 2);
M1_7_integrals = trapz(time, M1_7, 2); 
M1_8_integrals = trapz(time, M1_8, 2); 

scale1 = 1e-7; scale2 = 1e-9;
[f_M1_1, xi_M1_1] = ksdensity(M1_1_integrals.*scale1);
[f_M1_2, xi_M1_2] = ksdensity(M1_2_integrals.*scale1);
[f_M1_3, xi_M1_3] = ksdensity(M1_3_integrals.*scale1);
[f_M1_4, xi_M1_4] = ksdensity(M1_4_integrals.*scale1);
[f_M1_5, xi_M1_5] = ksdensity(M1_5_integrals.*scale1);
[f_M1_6, xi_M1_6] = ksdensity(M1_6_integrals.*scale1);
[f_M1_7, xi_M1_7] = ksdensity(M1_7_integrals.*scale1);
[f_M1_8, xi_M1_8] = ksdensity(M1_8_integrals.*scale1);

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
xi_data = {xi_M1_1, xi_M1_2, xi_M1_3, xi_M1_4, xi_M1_5, xi_M1_6, xi_M1_7, xi_M1_8};
f_data = {f_M1_1, f_M1_2, f_M1_3, f_M1_4, f_M1_5, f_M1_6, f_M1_7, f_M1_8};
for i = 1:length(xi_data)
    plot(xi_data{i}, f_data{i}, 'LineWidth', LineWidth, 'Color', Colors{i});
end
legend('Mouse 1', 'Mouse 2', 'Mouse 3', 'Mouse 4', 'Mouse 5', 'Mouse 6', 'Mouse 7', 'Mouse 8', ...
    'Location', 'northeast', 'NumColumns', 2, 'FontSize', FontSize-2);
xlabel('M1-type macrophages AUC (\times 10^7 cell)');
ylabel('Probability density');
title('Distribution curve');
set(gca, 'FontSize', FontSize, 'FontName', 'Arial');

subplot(1,3,2)
T_mean = [mean(T_1_integrals); mean(T_2_integrals); mean(T_3_integrals); ...
    mean(T_4_integrals); mean(T_5_integrals); mean(T_6_integrals); ...
    mean(T_7_integrals); mean(T_8_integrals)].*scale2;
M1_mean = [mean(M1_1_integrals); mean(M1_2_integrals); mean(M1_3_integrals); ...
    mean(M1_4_integrals); mean(M1_5_integrals); mean(M1_6_integrals); ...
    mean(M1_7_integrals); mean(M1_8_integrals)].*scale1;
hold on; grid on; box on;
for i = 1:8
    scatter(M1_mean(i), T_mean(i), MarkerSize+100, 'MarkerEdgeColor', [1 1 1], ...
        'MarkerFaceColor', Colors{i}, 'LineWidth', 1.5);
end
p = polyfit(M1_mean, T_mean, 1);
y_mean = polyval(p, M1_mean);
plot(M1_mean, y_mean, '-', 'Color', [0 0 0], 'LineWidth', 1.5);
[R_matrix, P_matrix] = corrcoef(M1_mean, T_mean);
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
xlim([2.75 5.25]); ylim([1 3]);
xlabel('M1-type macrophages AUC (\times 10^7 cell)');
ylabel('Tumor AUC (\times 10^9 mm^3)');
legend('Mouse 1', 'Mouse 2', 'Mouse 3', 'Mouse 4', 'Mouse 5', 'Mouse 6', 'Mouse 7', 'Mouse 8', ...
    'Location', 'northeast', 'NumColumns', 2, 'FontSize', FontSize-2);
set(gca, 'FontSize', FontSize, 'FontName', 'Arial');
title('Group');

subplot(1,3,3)
y_data = [T_1_integrals; T_2_integrals; T_3_integrals; T_4_integrals; ...
    T_5_integrals; T_6_integrals; T_7_integrals; T_8_integrals].*scale2;
x_data = [M1_1_integrals; M1_2_integrals; M1_3_integrals; M1_4_integrals; ...
    M1_5_integrals; M1_6_integrals; M1_7_integrals; M1_8_integrals].*scale1;
hold on; grid on; box on;
M1_data = {M1_1_integrals, M1_2_integrals, M1_3_integrals, M1_4_integrals, ...
    M1_5_integrals, M1_6_integrals, M1_7_integrals, M1_8_integrals};
T_data = {T_1_integrals, T_2_integrals, T_3_integrals, T_4_integrals, ...
          T_5_integrals, T_6_integrals, T_7_integrals, T_8_integrals};
for i = 1:length(M1_data)
    x = M1_data{i} .* scale1;
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
xlabel('M1-type macrophages AUC (\times 10^7 cell)');
ylabel('Tumor AUC (\times 10^9 mm^3)');
title('Individual');
xlim([1.5 13.5]); ylim([1.2 3.5]);
set(gca, 'FontSize', FontSize, 'FontName', 'Arial');

set(gcf,'unit','centimeters','position',[3 10 45 10]);
print('Figure/Fig8ABC','-dpng','-r200');
