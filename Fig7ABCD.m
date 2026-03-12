function Fig7ABCD()

% close all; clc; clear

% Load data
Q = 2501;
T = cell2mat(struct2cell(load('./Data/Tumor.mat')));
T_1 = T(:,1:Q,1); T_2 = T(:,1:Q,2); T_3 = T(:,1:Q,3); T_4 = T(:,1:Q,4);
T_5 = T(:,1:Q,5); T_6 = T(:,1:Q,6); T_7 = T(:,1:Q,7); T_8 = T(:,1:Q,8);
M1 = cell2mat(struct2cell(load('./Data/M1.mat')));
M1_1 = M1(:,1:Q,1); M1_2 = M1(:,1:Q,2); M1_3 = M1(:,1:Q,3); M1_4 = M1(:,1:Q,4);
M1_5 = M1(:,1:Q,5); M1_6 = M1(:,1:Q,6); M1_7 = M1(:,1:Q,7); M1_8 = M1(:,1:Q,8);
M2 = cell2mat(struct2cell(load('./Data/M2.mat')));
M2_1 = M2(:,1:Q,1); M2_2 = M2(:,1:Q,2); M2_3 = M2(:,1:Q,3); M2_4 = M2(:,1:Q,4);
M2_5 = M2(:,1:Q,5); M2_6 = M2(:,1:Q,6); M2_7 = M2(:,1:Q,7); M2_8 = M2(:,1:Q,8);
M3 = cell2mat(struct2cell(load('./Data/M3.mat')));
M3_1 = M3(:,1:Q,1); M3_2 = M3(:,1:Q,2); M3_3 = M3(:,1:Q,3); M3_4 = M3(:,1:Q,4);
M3_5 = M3(:,1:Q,5); M3_6 = M3(:,1:Q,6); M3_7 = M3(:,1:Q,7); M3_8 = M3(:,1:Q,8);

% Data processing
time = 0:0.01:25;
T_1_integrals = trapz(time, T_1, 2); std_value_T_1 = std(T_1_integrals);
T_2_integrals = trapz(time, T_2, 2); std_value_T_2 = std(T_2_integrals);
T_3_integrals = trapz(time, T_3, 2); std_value_T_3 = std(T_3_integrals);
T_4_integrals = trapz(time, T_4, 2); std_value_T_4 = std(T_4_integrals);
T_5_integrals = trapz(time, T_5, 2); std_value_T_5 = std(T_5_integrals);
T_6_integrals = trapz(time, T_6, 2); std_value_T_6 = std(T_6_integrals);
T_7_integrals = trapz(time, T_7, 2); std_value_T_7 = std(T_7_integrals);
T_8_integrals = trapz(time, T_8, 2); std_value_T_8 = std(T_8_integrals);

integral = [mean(T_1_integrals); mean(T_2_integrals); mean(T_3_integrals);
    mean(T_4_integrals); mean(T_5_integrals); mean(T_6_integrals); 
    mean(T_7_integrals); mean(T_8_integrals)];
std_value = [std_value_T_1; std_value_T_2; std_value_T_3;
    std_value_T_4; std_value_T_5; std_value_T_6;
    std_value_T_7; std_value_T_8];

name = {'Mouse 1', 'Mouse 2', 'Mouse 3', 'Mouse 4', 'Mouse 5', 'Mouse 6', 'Mouse 7', 'Mouse 8'};
M1_data = {M1_1, M1_2, M1_3, M1_4, M1_5, M1_6, M1_7, M1_8};
M2_data = {M2_1, M2_2, M2_3, M2_4, M2_5, M2_6, M2_7, M2_8};
M3_data = {M3_1, M3_2, M3_3, M3_4, M3_5, M3_6, M3_7, M3_8};

Color = ColorMatrix();
FontSize = 14;
MarkerSize = 50;
LineWidth = 2;
Colors = {cell2mat(Color(3,9)); cell2mat(Color(3,10)); cell2mat(Color(3,1)); ...
    cell2mat(Color(5,1)); cell2mat(Color(7,1)); cell2mat(Color(3,11)); ...
    cell2mat(Color(3,13)); cell2mat(Color(5,13))};

figure()

subplot(2,2,1)
hold on; grid on; box on;
h = bar(name, integral*1e-9);
h.FaceColor = 'flat';
h.CData = [cell2mat(Color(3,9)); cell2mat(Color(3,10)); cell2mat(Color(3,1)); ...
    cell2mat(Color(5,1)); cell2mat(Color(7,1)); cell2mat(Color(3,11)); ...
    cell2mat(Color(3,13)); cell2mat(Color(5,13))];
errorbar(1:8, integral*1e-9, std_value*1e-9, 'k.', 'LineWidth', LineWidth-1, 'CapSize', 15);
yline(2.4, 'r--', 'LineWidth', LineWidth, 'FontSize', 10);
yline(1.8, 'b--', 'LineWidth', LineWidth, 'FontSize', 10);
ylim([1 3.5]);
ylabel('Tumor AUC (\times 10^9 mm^3)');
set(gca, 'FontSize', FontSize, 'FontName', 'Arial');

subplot(2,2,2)
hold on; grid on; box on;
for i = 1:length(M1_data)
    plot(time, mean(M1_data{i}, 1), '-o', ...
        'Color', Colors{i}, 'LineWidth', LineWidth, ...
        'MarkerIndices', 1:300:length(time), ... 
        'MarkerSize', 6, ...
        'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', Colors{i});
end
ylabel('Mean M1-type dynamics');
set(gca, 'FontSize', FontSize, 'FontName', 'Arial');
legend('Mouse 1', 'Mouse 2', 'Mouse 3', 'Mouse 4', 'Mouse 5', 'Mouse 6', 'Mouse 7', 'Mouse 8', ...
    'Location', 'northeast', 'NumColumns', 2);
xlabel('Time (days)');
xlim([0 25]); ylim([0 4e6]);

subplot(2,2,3)
hold on; grid on; box on;
for i = 1:length(M3_data)
    plot(time, mean(M3_data{i}, 1), '-s', ...
        'Color', Colors{i}, 'LineWidth', LineWidth, ...
        'MarkerIndices', 1:300:length(time), ... 
        'MarkerSize', 6, ...
        'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', Colors{i});
end
ylabel('Mean M3-type dynamics');
set(gca, 'FontSize', FontSize, 'FontName', 'Arial');
legend('Mouse 1', 'Mouse 2', 'Mouse 3', 'Mouse 4', 'Mouse 5', 'Mouse 6', 'Mouse 7', 'Mouse 8', ...
    'Location', 'northeast', 'NumColumns', 2);
xlabel('Time (days)');
xlim([0 25]); ylim([0 2.5e6]);

subplot(2,2,4)
hold on; grid on; box on;
for i = 1:length(M2_data)
    plot(time, mean(M2_data{i}, 1), '-^', ...
        'Color', Colors{i}, 'LineWidth', LineWidth, ...
        'MarkerIndices', 1:300:length(time), ... 
        'MarkerSize', 6, ...
        'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', Colors{i});
end
ylabel('Mean M2-type dynamics');
set(gca, 'FontSize', FontSize, 'FontName', 'Arial');
legend('Mouse 1', 'Mouse 2', 'Mouse 3', 'Mouse 4', 'Mouse 5', 'Mouse 6', 'Mouse 7', 'Mouse 8', ...
    'Location', 'southeast', 'NumColumns', 2);
xlabel('Time (days)');
xlim([0 25]); ylim([0 4e6]);

set(gcf,'unit','centimeters','position',[3 10 30 24])
print('Figure/Fig7ABCD','-dpng','-r200')