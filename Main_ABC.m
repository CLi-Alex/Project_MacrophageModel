function Main_ABC()

% Input model parameters
RR = 0.5;

% Load data
Data = [8 11 13 15 18 20 23; 86 268.02 430.35 670.93 1103.47 1556.73 1732.68;
    8 11 13 15 18 20 23; 74.66 284.71 430.71 619.76 1033.08 1350.51 1705.42;
    8 11 13 15 18 20 23; 72 200.29 326.36 610.41 1267.32 1646 2645.34;
    8 11 13 15 18 20 23; 71.54 195.76 351.1 680.19 1171.53 1716.05 2412.35;
    8 11 13 15 18 20 NaN; 62 211.76 454.96 840.17 1960.9 2388.14 NaN;
    8 11 13 15 18 20 23; 47.7 190.33 371.4 608.32 860.8 1265.47 1807.26;
    8 11 13 15 18 20 23; 39.64 145.12 214.51 340.28 540.06 893.15 1340.64;
    8 11 13 15 18 20 23; 23.88 86.49 164.9 315 596.66 1015.85 1705.9];

% Parameter sampling
m = 1500;
mm = 3101;
Num = 100;
Matrix = zeros(8,Num,8);
R = zeros(8,Num);
W = zeros(8,Num);
Tumor_prior = zeros(m, mm, 8);  Tumor = zeros(Num, mm, 8);
M1_prior    = zeros(m, mm, 8);  M1    = zeros(Num, mm, 8);
M2_prior    = zeros(m, mm, 8);  M2    = zeros(Num, mm, 8);
M3_prior    = zeros(m, mm, 8);  M3    = zeros(Num, mm, 8);
for k = 1:8
    Matrix_ABC = ABC_Paras(k,m);
    % Approximate Bayesian computation
    W_ABC = zeros(1,m);
    R_ABC = zeros(1,m);
    for j = 1:m
        % Parameter assignment
        M = Matrix_ABC(:,j);
        Y = ABC_Model(M);
        if k == 5
            D = Data(2*k-1:2*k,1:6);
        else
            D = Data(2*k-1:2*k,1:7);
        end
        T = 1 + 100 .* D(1,:);
        P = D(2,:);
        Q = Y(1,T).*6.95e-6;
        R_ABC(j) = Method_R2(P, Q);
        if R_ABC(j) > RR
            Indicator = 1;
        else
            Indicator = 0;
        end
        W_ABC(1,j) = Indicator;
        Tumor_prior(j,:,k) = Y(1,:);
        M1_prior(j,:,k)    = Y(2,:);
        M2_prior(j,:,k)    = Y(3,:);
        M3_prior(j,:,k)    = Y(4,:);
    end
    [~, topNum_idx] = maxk(R_ABC, Num);
    Matrix(:,:,k) = Matrix_ABC(:,topNum_idx);
    R(k,:) = R_ABC(:,topNum_idx);
    W(k,:) = W_ABC(:,topNum_idx);
    Tumor(:,:,k) = Tumor_prior(topNum_idx, :,k);
    M1(:,:,k)    = M1_prior(topNum_idx, :,k);
    M2(:,:,k)    = M2_prior(topNum_idx, :,k);
    M3(:,:,k)    = M3_prior(topNum_idx, :,k);
end

save('Data/Matrix.mat', 'Matrix');
save('Data/R.mat', 'R');
save('Data/W.mat', 'W');
save('Data/Tumor.mat', 'Tumor');
save('Data/M1.mat', 'M1');
save('Data/M2.mat', 'M2');
save('Data/M3.mat', 'M3');