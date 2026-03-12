function main()
    % Main function for tumor-macrophage interaction model
    % Solves ODE system using Euler method
    
    % Load parameters
    params = parameters();
    
    % Set up time vector
    tspan = params.t_start:params.dt:params.t_end;
    n_steps = length(tspan);
    
    % Initialize storage arrays
    Y = zeros(4, n_steps);
    Y(:, 1) = params.y0';
    
    % Euler method integration
    for i = 1:n_steps-1
        Y(:, i+1) = Y(:, i) + params.dt * ode_system(tspan(i), Y(:, i), params);
        
        % Ensure non-negative values
        % Y(:, i+1) = max(Y(:, i+1), 0);
    end
    
    % Extract results
    C = Y(1, :);
    M1 = Y(2, :);
    M2 = Y(3, :);
    M3 = Y(4, :);
    
    % Plot results
    Fig2(tspan, C, M1, M2, M3);
end

