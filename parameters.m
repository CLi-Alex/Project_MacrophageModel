function params = parameters()
    % Parameters for tumor-macrophage interaction model
    % Returns parameter structure
    
    close all;

    % Tumor growth parameters
    params.alpha = 0.431;     % Tumor growth rate
    params.delta = 1e-9;      % Promotion effect of M2/M3 on tumor growth
    params.K = 3e8;           % Tumor carrying capacity (3e8 Lai.BullMathBiol.2025)
    params.eta = 1.1e-8;      % Killing rate of M1/M3 on tumor
    params.p = 0.5;           % Proportion of M3 that promotes tumor growth
    
    % Macrophage transformation parameters
    params.kappa = 0.4;      % M0 polarization rate
    params.beta13 = 0.5;     % M1->M3 conversion rate
    params.beta31 = 0;       % M3->M1 conversion rate
    params.beta32 = 0.5;     % M3->M2 conversion rate
    params.beta23 = 0;       % M2->M3 conversion rate
    
    % Macrophage death parameters
    params.d1 = 0.056;      % M1 death rate
    params.d2 = 0.056;      % M2 death rate
    params.d3 = 0.056;      % M3 death rate
    
    % M2 polarization parameters
    params.gamma  = 0.1;    % M0->M2 basal polarization rate
    params.lambda = 0.5;    % Tumor promotion effect on M2 polarization
    
    % Hill function parameters
    params.K0 = 3e7;       % M2 polarization half-saturation constant
    params.K1 = 3e6;       % M1->M3 conversion half-saturation constant
    params.K2 = 3e6;       % M3->M2 conversion half-saturation constant
    
    % Constant population
    params.M0 = 4e5;       % Constant M0 population
    
    % Simulation parameters
    params.t_start = 0;     % Start time
    params.t_end   = 31;    % End time
    params.dt      = 0.01;  % Time step
    
    % Initial conditions [C, M1, M2, M3]
    params.y0 = [2e5, 4e6, 0, 0];  

end