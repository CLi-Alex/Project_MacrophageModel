function Fig2(t, C, M1, M2, M3)
    % Plot simulation results
    
    figure('Position', [100, 100, 1200, 800]);
    Color = ColorMatrix();

    x = [8 11 13 15 18 20 23];
    y = [59.6775 197.81 343.03625 585.6325 1066.7275 1478.9875 1907.08429];
    V = 6.95e-9 .* C;
    R2 = Method_R2(y, V(x.*100)*1e3);

    MR1 = M1 ./ (M1 + M2 + M3);
    MR2 = M2 ./ (M1 + M2 + M3);
    MR3 = M3 ./ (M1 + M2 + M3);
    
    LineWidth = 2;
    MarkerSize = 10;
    % Individual subplots
    subplot(1, 3, 1);
    hold on; grid on; box on;
    plot(t, V, 'r-', 'LineWidth', 2.5);
    scatter(x, y*1e-3, 85, 'o', 'MarkerFaceColor','r')
    text(0.1, 0.6, ['R^2=', num2str(R2,'%.3g')], ...
    'FontSize',16,'FontWeight','bold', 'Units', 'normalized')
    xlabel('Time (days)');
    ylabel('Tumor volume (\times 10^3 mm^3)');
    legend('Model simulation', 'Mean experimental data', 'Location', 'northwest','FontSize',14)
    set(gca,'FontSize',12);
    set(gca,'FontName','Arial');
    xlim([0 25]);

    subplot(1, 3, 2);
    hold on; grid on; box on;
    plot(t, M1, '-o', 'Color', 'b', 'LineWidth', LineWidth, 'MarkerIndices', 1:250:length(t), ... 
        'MarkerSize', MarkerSize, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'w');
    plot(t, M2, '-^', 'Color', 'r', 'LineWidth', LineWidth, 'MarkerIndices', 1:250:length(t), ... 
        'MarkerSize', MarkerSize, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'w');
    plot(t, M3, '-s', 'Color', cell2mat(Color(3,13)), 'LineWidth', LineWidth, 'MarkerIndices', 1:200:length(t), ... 
        'MarkerSize', MarkerSize, 'MarkerFaceColor', cell2mat(Color(3,13)), 'MarkerEdgeColor', 'w');
    xlabel('Time (days)');
    ylabel('Macrophages');
    legend('M1-type macrophages', 'M2-type macrophages', 'M3-type macrophages', 'Location', 'northeast','FontSize',12)
    set(gca,'FontSize',12);
    set(gca,'FontName','Arial');
    ylim([0 5e6]);
    xlim([0 25]);

    subplot(1, 3, 3);
    hold on; grid on; box on;
    plot(t, MR1, '-o', 'Color', 'b', 'LineWidth', LineWidth, 'MarkerIndices', 1:250:length(t), ... 
        'MarkerSize', MarkerSize, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'w');
    plot(t, MR2, '-^', 'Color', 'r', 'LineWidth', LineWidth, 'MarkerIndices', 1:250:length(t), ... 
        'MarkerSize', MarkerSize, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'w');
    plot(t, MR3, '-s', 'Color', cell2mat(Color(3,13)), 'LineWidth', LineWidth, 'MarkerIndices', 1:200:length(t), ... 
        'MarkerSize', MarkerSize, 'MarkerFaceColor', cell2mat(Color(3,13)), 'MarkerEdgeColor', 'w');
    xlabel('Time (days)');
    ylabel('% Macrophages');
    set(gca,'FontSize',12);
    set(gca,'FontName','Arial');
    xlim([0 25]);

    % Image output Settings
    set(gcf,'unit','centimeters','position',[3 10 45 8]);
    print('-dpng','-r200','Figure/Fig2.png');
end