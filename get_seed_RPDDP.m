%%=========================================================================
% Brief: Get the initial seeds using PCA
%   Input:
%       Covariance Image
%   Output
%       seed, data1, data2
% Autor: Naiallen Carvalho
% Computação Aplicada
% Instituto Nacional de Pesquisas Espaciais - INPE
% 2020
%==========================================================================

function [ seed, data_cov1, data_cov2] = get_seed_RPDDP(imacov_vec)

%% Compute parameters -------------------------------------------------
% data = [sqrt(imacov_vec(:,1)) sqrt(imacov_vec(:,5)) sqrt(imacov_vec(:,9))];
% dataz = abs(zscore(data));
% index = zeros(size(data,1),1);
% index( (dataz(:,1)>2.5) | (dataz(:,2)>2.5) | (dataz(:,3)>2.5),:) = 1;
% imacov_vec(index==1,:) = [];

[~, n_bands] = size(imacov_vec);
Sigma = intrisic_mean( imacov_vec );%average_covariance(imacov_vec);
m = diag(Sigma)';

%% Compue eigenvalue and eigenvector---------------------------------
[V, D] = eig(Sigma);

%% Sort eigenvalue---------------------------------------------------
[eigenval, b] = sort(diag(D), 'descend');

%% Sort eigenvector--------------------------------------------------
eigenvec = zeros(size(V));
for ev = 1:size(eigenval,1)
    eigenvec(:, ev) = V(:,b(ev));
end

%% Find the signal
d = zeros(size(imacov_vec,1),1);
ab = zeros(size(imacov_vec,1),1);
if (n_bands == 9)
    data = [imacov_vec(:, 1) imacov_vec(:, 5) imacov_vec(:, 9)];
elseif(n_bands == 18)
    data = [imacov_vec(:, 1) imacov_vec(:, 5) imacov_vec(:, 9) imacov_vec(:, 10) imacov_vec(:, 14) imacov_vec(:, 18)];
end
for ii = 1:size(imacov_vec,1)   
    ab(ii) = dot(eigenvec(:, 1), data(ii,:)-m);
    d(ii) = acosd( dot(eigenvec(:, 1), data(ii,:)-m)/ ( norm(eigenvec(:, 1))*norm(data(ii,:)-m)) );
end 
d(d>90) = -1;

data_cov1 = imacov_vec(d>0,:);
data_cov2 = imacov_vec(d<0,:);
% 
% %% Compute the normal plane to the biggest eigenvalue----------------
% plan = cross(eigenvec(:, 2), eigenvec(:, 3));
% 
% %% Find plan equation coefficients-----------------------------------
% A = plan(1);
% B = plan(2);
% C = plan(3);
% D = -dot(plan,m);
% 
% %% Divide data given the plan----------------------------------------------
% d = A*sqrt(imacov_vec(:, 1)) + B*sqrt(imacov_vec(:, 5)) + C*sqrt(imacov_vec(:, 9)) + D;
% data_cov21 = imacov_vec(d>0,:);
% data_cov22 = imacov_vec(d<0,:);

%% Get Seed 1--------------------------------------------------------------
m_cov = intrisic_mean(  imacov_vec(d>0,:) );
seed{1} = m_cov;

%% Get Seed 2--------------------------------------------------------------
m_cov = intrisic_mean( imacov_vec(d<0,:) );
seed{2} =  m_cov;

end

