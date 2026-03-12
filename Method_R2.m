%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is used to calculate R^2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function R2 = Method_R2(observed, simulated)

observed = observed(:);
simulated = simulated(:);

RSS = sum((observed - simulated).^2);      
TSS = sum((observed - mean(observed)).^2); 

R2 = 1 - (RSS / TSS);

end
