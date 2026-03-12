function Fig3()

clear; clc;

param_names = {'\alpha','\it K','\delta','\it p','\eta','\kappa','\gamma',...
    '\lambda','\it M_0','\beta_{13}','\beta_{32}','\it K_0','\it K_1',...
    '\it K_2','\it d_1','\it d_{2}','\it d_3'};

load('Data/Si_C.mat'); load('Data/STi_C.mat'); load('Data/Si_M1.mat'); load('Data/STi_M1.mat');
load('Data/Si_M2.mat'); load('Data/STi_M2.mat'); load('Data/Si_M3.mat'); load('Data/STi_M3.mat');

time = 1:25;
Color = ColorMatrix();
colors = [cell2mat(Color(2,1)); cell2mat(Color(2,4)); cell2mat(Color(2,6));
    cell2mat(Color(2,8)); cell2mat(Color(2,9)); cell2mat(Color(2,10)); 
    cell2mat(Color(2,11)); cell2mat(Color(2,13)); cell2mat(Color(2,14));
    cell2mat(Color(6,1)); cell2mat(Color(6,4)); cell2mat(Color(6,6)); 
    cell2mat(Color(6,8)); cell2mat(Color(6,9)); cell2mat(Color(6,10)); 
    cell2mat(Color(6,11)); cell2mat(Color(6,13)); cell2mat(Color(6,14))
];

markers = {'o', 'o', 'o', 's', 's', 's', 'd', 'd', 'd', '^', '^', '^', 'p', 'p', 'p', 'v', 'v'};
line_styles = {'-', '--', ':', '-.', '-', '--', ':', '-.', '-', '--', ':', '-.', '-', '--', ':', '-.', '-'};

figure();
tlayout = tiledlayout(1, 4, 'TileSpacing', 'loose', 'Padding', 'loose');
subplot_titles = {'Tumor cells', 'M1-type macrophages', 'M2-type macrophages', 'M3-type macrophages'};
data_to_plot = {STi_C, STi_M1, STi_M2, STi_M3};

plot_handles = []; 

for plot_idx = 1:4
    nexttile;
    hold on; grid on;
    for i = 1:17
        h = plot(time, data_to_plot{plot_idx}(i,:), ...
            'Color', colors(i,:), ...
            'LineStyle', line_styles{i}, ...
            'LineWidth', 1.8, ...
            'Marker', markers{i}, ...
            'MarkerSize', 3.5, ...
            'MarkerEdgeColor', colors(i,:), ...
            'MarkerFaceColor', colors(i,:));
    end
    
    hold off;
    set(gca, 'FontSize', 14, 'LineWidth', 1, 'TickDir', 'in', 'Box', 'on');
    xlabel('Time (days)', 'FontSize', 14);
    ylabel('Total order index', 'FontSize', 14);
    title(subplot_titles{plot_idx}, 'FontSize', 14);
    xlim([0, 25]);
    if plot_idx == 1
        legend('\alpha','\it K','\delta','\it p','\eta','\kappa','\gamma', ...
    '\lambda','\it M_0','\beta_{13}','\beta_{32}','\it K_0','\it K_1', ...
    '\it K_2','\it d_1','\it d_{2}','\it d_3', ...
       'Interpreter', 'tex','FontSize', 12, 'NumColumns', 3, ...           
       'Orientation', 'horizontal', 'Box', 'off', 'Location', 'best');
    end
end

set(gcf,'unit','centimeters','position',[3 10 45 9]);
print('Figure/Fig3','-dpng','-r200');