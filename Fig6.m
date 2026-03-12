function Fig6()

% close all; clear; clc

% Load data
Tumor = cell2mat(struct2cell(load('./Data/Tumor.mat')));
R = cell2mat(struct2cell(load('./Data/R.mat')));

% Experimental data
ExpData = [8 11 13 15 18 20 23; 86 268.02 430.35 670.93 1103.47 1556.73 1732.68;
    8 11 13 15 18 20 23; 74.66 284.71 430.71 619.76 1033.08 1350.51 1705.42;
    8 11 13 15 18 20 23; 72 200.29 326.36 610.41 1267.32 1646 2645.34;
    8 11 13 15 18 20 23; 71.54 195.76 351.1 680.19 1171.53 1716.05 2412.35;
    8 11 13 15 18 20 NaN; 62 211.76 454.96 840.17 1960.9 2388.14 NaN;
    8 11 13 15 18 20 23; 47.7 190.33 371.4 608.32 860.8 1265.47 1807.26;
    8 11 13 15 18 20 23; 39.64 145.12 214.51 340.28 540.06 893.15 1340.64;
    8 11 13 15 18 20 23; 23.88 86.49 164.9 315 596.66 1015.85 1705.9];

% Figure 1
figure()
Color = ColorMatrix();
time = 0:0.01:31;
LineWidth = 1;
FontSize = 12;
MarkerSize = 50;
nbins = 15;
ColorsLine = {cell2mat(Color(1,9)); cell2mat(Color(1,10)); cell2mat(Color(1,1)); ...
    cell2mat(Color(3,1)); cell2mat(Color(5,1)); cell2mat(Color(1,11)); ...
    cell2mat(Color(1,13)); cell2mat(Color(5,13))};
ColorsPoint = {cell2mat(Color(3,9)); cell2mat(Color(3,10)); cell2mat(Color(3,1)); ...
    cell2mat(Color(5,1)); cell2mat(Color(7,1)); cell2mat(Color(3,11)); ...
    cell2mat(Color(1,13)); cell2mat(Color(5,13))};

for i = 1:8
    subplot(2,4,i)
    hold on; box on; grid on;
    [~, pos] = max(R(i,:));
    Data = Tumor(:,:,i).*6.95e-9;
    Max = max(Data, [], 1); Min = min (Data, [], 1);
    plot(time, Data(pos,:), '-', 'LineWidth', LineWidth+1, 'Color', ColorsLine{i});
    fill([time, fliplr(time)], [Max, fliplr(Min)], ...
        ColorsPoint{i}, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    if i == 5
        scatter(ExpData(2*i-1,1:6), ExpData(2*i,1:6)*1e-3, MarkerSize, ...
            'MarkerEdgeColor', ColorsLine{i}, 'MarkerFaceColor', ColorsPoint{i}, ...
            'LineWidth', LineWidth+1);
    else
        scatter(ExpData(2*i-1,:), ExpData(2*i,:)*1e-3, MarkerSize, ...
            'MarkerEdgeColor', ColorsLine{i}, 'MarkerFaceColor', ColorsPoint{i}, ...
            'LineWidth', LineWidth+1);
    end
    plot(time, Max, ':', 'LineWidth', LineWidth, 'Color', ColorsLine{i});
    plot(time, Min, ':', 'LineWidth', LineWidth, 'Color', ColorsLine{i});
    text(0.1, 0.5, sprintf('R^2 = %.3f', R(i,pos)), 'Units', 'normalized', ...
        'FontSize', FontSize+2, 'FontWeight', 'bold');
    xlim([0,25]);
    ylim([0,3.25]);
    xlabel('Time (days)');
    ylabel('Tumor volume (\times 10^3 mm^3)');
    title(sprintf('Mouse %d', i))
    set(gca,'FontSize', FontSize);
    set(gca,'FontName','Arial');
    legend('Best fit', 'Posterior inference', 'Data','Location','northwest');
end

set(gcf,'unit','centimeters','position',[3 10 48 18])
print('Figure/Fig6A','-dpng','-r200')

% Figure 2
figure()
for i = 1:8
    subplot(2,4,i)
    hold on; box on; grid on;
    histogram(R(i,:),nbins, 'FaceColor', ColorsPoint{i})
    xlabel('R^2');
    ylabel('Count');
    title(sprintf('Mouse %d', i))
    set(gca,'FontSize', FontSize);
    set(gca,'FontName','Arial');
end

set(gcf,'unit','centimeters','position',[3 10 48 10])
print('Figure/Fig6B','-dpng','-r200')