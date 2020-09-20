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
function [ seed, data1, data2 ] = get_seeds_random( imacov_vec )

s = size(imacov_vec, 1);

index1 = randi(s);
index2 = randi(s);

data1 = imacov_vec(index1,:);
data2 = imacov_vec(index2,:);


%% Get Seed 1--------------------------------------------------------------
seed{1} = reshape(imacov_vec(index1,:), 3, 3)';

%% Get Seed 2--------------------------------------------------------------
seed{2} =  reshape(imacov_vec(index2,:), 3, 3)';

end

