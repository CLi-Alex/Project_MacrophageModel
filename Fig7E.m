function Fig7E()

% Load data
Matrix = cell2mat(struct2cell(load('./Data/Matrix.mat')));

group_def = {
    [3, 4, 5],   
    [1, 2, 6],  
    [7, 8]    
};
n_groups = length(group_def);
param_names = {'alpha', 'K', 'M0', 'Kappa', 'Beta13', 'Beta32', 'Eta', 'P'};
n_params = length(param_names);

means_all = zeros(n_params, n_groups);
stds_all = zeros(n_params, n_groups);

for param_idx = 1:n_params
    for group_idx = 1:n_groups
        m_indices = group_def{group_idx};
        group_data = [];
        for m_idx = m_indices
            current_data = Matrix(param_idx, :, m_idx);
            group_data = [group_data, current_data];
        end
        means_all(param_idx, group_idx) = mean(group_data);
        stds_all(param_idx, group_idx) = std(group_data);
    end
end

figure('Position', [100, 100, 1600, 800]);
Color = ColorMatrix();
FontSize = 16;
bar_colors = [cell2mat(Color(1,1)); cell2mat(Color(1,11)); cell2mat(Color(1,13))];

ylimits = [
    0.37 0.47;      % alpha
    2e8 5e8;        % K
    0 1e6;          % M0
    0 1;            % Kappa
    0 1;            % Beta13
    0 1;            % Beta32
    0 3.6e-8;       % Eta
    0 1             % P
];

titles = {'\alpha', 'K', 'M_0', '\kappa','\beta_{13}', '\beta_{32}', '\eta', 'p'};

for i = 1:n_params
    subplot(4, 2, i);
    hold on; grid on; box on;
    x_categories = {'I', 'II', 'III'};
    h = bar(x_categories, means_all(i, :), ...
        'FaceColor', 'none', ...
        'LineWidth', 2.5, ...
        'EdgeColor', 'flat');
    h.CData = bar_colors;
    errorbar(1:3, means_all(i, :), stds_all(i, :), ...
        'k.', 'LineWidth', 1.5, 'CapSize', 15);
    set(gca, 'FontSize', FontSize-4, 'FontName', 'Arial');
    if ~any(isnan(ylimits(i, :)))
        ylim(ylimits(i, :));
        num_ticks = 3;
        yticks(linspace(ylimits(i, 1), ylimits(i, 2), num_ticks));
    end
    title(titles{i}, 'FontSize', FontSize, 'Interpreter', 'tex');
end

set(gcf, 'unit', 'centimeters', 'position', [3 10 18 26]);
print('Figure/Fig7E','-dpng','-r200')
