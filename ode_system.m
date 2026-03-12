function dydt = ode_system(t, y, params)
    % ODE system definition
    % y = [C; M1; M2; M3]
    
    C = y(1);
    M1 = y(2);
    M2 = y(3);
    M3 = y(4);
    
    % Hill functions
    hill_M2_polarization = C / (params.K0 + C);
    hill_M1_to_M3 = C / (params.K1 + C);
    hill_M3_to_M2 = C / (params.K2 + C);
    
    % Tumor cell equation
    dC_dt = params.alpha * (1 + params.delta * (M2 + params.p * M3)) * ...
            (1 - C/params.K) * C - params.eta * (M1 + (1 - params.p) * M3) * C;
    
    % M1 macrophage equation
    dM1_dt = params.kappa * params.M0 - params.beta13 * hill_M1_to_M3 * M1 + ...
             params.beta31 * M3 - params.d1 * M1;
    
    % M2 macrophage equation
    dM2_dt = params.gamma * (1 + params.lambda * hill_M2_polarization) * params.M0 + ...
             params.beta32 * hill_M3_to_M2 * M3 - params.beta23 * M2 - params.d2 * M2;
    
    % M3 macrophage equation
    dM3_dt = params.beta13 * hill_M1_to_M3 * M1 - params.beta31 * M3 - ...
             params.beta32 * hill_M3_to_M2 * M3 + params.beta23 * M2 - params.d3 * M3;
    
    dydt = [dC_dt; dM1_dt; dM2_dt; dM3_dt];
end

