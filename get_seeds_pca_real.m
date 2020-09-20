%%=========================================================================
% Brief: Get the initial seeds using PCA
%   Input:
%       Covariance Image: 9 bands image
%   Output
%       seed, data1, data2
% Autor: Naiallen Carvalho
% Computação Aplicada
% Instituto Nacional de Pesquisas Espaciais - INPE
% 2020
%==========================================================================

%%
function [ seed, data1, data2 ] = get_seeds_pca_real( data )


%% Compute parameters -------------------------------------------------
m = mean(data);
Sigma = get_covariance(data);

%% Compute eigenvalue and eigenvector---------------------------------
[V, D] = eig(Sigma);

%% Sort eigenvalue---------------------------------------------------
[eigenval, b] = sort(diag(D), 'descend');

%% Sort eigenvector--------------------------------------------------
eigenvec = zeros(size(V));
for ev = 1:size(eigenval,1)
    eigenvec(:, ev) = V(:,b(ev));
end

%% Divide data given the plan----------------------------------------------
d = zeros(size(data,1),1);
ab = zeros(size(data,1),1);
for ii = 1:size(data,1)   
    ab(ii) = dot(eigenvec(:, 1), data(ii,:)-m);
    d(ii) = acosd( dot(eigenvec(:, 1), data(ii,:)-m)/ ( norm(eigenvec(:, 1))*norm(data(ii,:)-m)) );
end 
d(d>90) = -1;

data1 = data(d>0,:);
data2 = data(d<0,:);
%% Get Seed 1--------------------------------------------------------------
seed{1}.cov = get_covariance(data1);
seed{1}.m = mean(data1);

%% Get Seed 2--------------------------------------------------------------
seed{2}.cov = get_covariance(data2);
seed{2}.m = mean(data2);
end

