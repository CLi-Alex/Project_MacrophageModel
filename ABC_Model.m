function Y = ABC_Model(M)

% Load parameters
params = parameters();
params.alpha  = M(1); 
params.K      = M(2);  
params.M0     = M(3);
params.kappa  = M(4);
params.beta13 = M(5);
params.beta32 = M(6);
params.eta    = M(7);   
params.p      = M(8);

% Set up time vector
tspan = params.t_start:params.dt:params.t_end;
n_steps = length(tspan);

% Initialize storage arrays
Y = zeros(4, n_steps);
Y(:, 1) = params.y0';

% Euler method integration
for i = 1:n_steps-1
    Y(:, i+1) = Y(:, i) + params.dt * ode_system(tspan(i), Y(:, i), params);
end

end