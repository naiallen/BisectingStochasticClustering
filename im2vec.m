%%=========================================================================
% Brief: Convert a image to a vector form
%   Input:
%       image: [nrow, ncol, n_band]
%   Output
%       vector: [nrow*ncol, n_band]
% Autor: Naiallen Carvalho
% Computação Aplicada
% Instituto Nacional de Pesquisas Espaciais - INPE
% 2020
%==========================================================================
function [ vector ] = im2vec( image)
[~, ~, n_bands] = size(image);
vector = [];
for ii =1:n_bands
    clear aux
    aux =  image(:, :, ii); aux=(aux(~isnan(aux)));
    vector = [vector aux(:)];
end
end

