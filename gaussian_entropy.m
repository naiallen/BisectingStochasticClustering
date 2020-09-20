%% ========================================================================
% Brief: Compute the Multivariate Normal Entropy
%   Input:
%       parameters
%   Output
%       H: entropy
% Autor: Naiallen Carvalho
% Computação Aplicada
% Instituto Nacional de Pesquisas Espaciais - INPE
% 2020
%==========================================================================

function H = gaussian_entropy( parameters)
N = parameters.q;
H = 0.5 * log( ((2*pi*exp(1))^N)*det(parameters.cov) );
end

