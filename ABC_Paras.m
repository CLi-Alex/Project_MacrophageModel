function Matrix = ABC_Paras(i,m)

Matrix = zeros(8,m);
rng(i);

Matrix(1,:) = 0.38   + (0.46 - 0.38)  * rand(1,m); % Alpha
Matrix(2,:) = 2.17e8 + (5e8 - 2.17e8) * rand(1,m); % K
Matrix(3,:) = 1e5    + (1e6 - 1e5)    * rand(1,m); % M0
Matrix(4,:) = 0      + (1 - 0)        * rand(1,m); % Kappa
Matrix(5,:) = 0      + (1 - 0)        * rand(1,m); % Beta_13
Matrix(6,:) = 0      + (1 - 0)        * rand(1,m); % Beta_32
Matrix(7,:) = 1e-8   + (3e-8 - 1e-8)  * rand(1,m); % Eta
Matrix(8,:) = 0      + (1 - 0)        * rand(1,m); % P

end