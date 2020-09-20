%%=========================================================================
% Brief: Compute the Entropy
%   Input:
%       parameters
%   Output
%       H: entropy
% Autor: Naiallen Carvalho
% Computação Aplicada
% Instituto Nacional de Pesquisas Espaciais - INPE
% 2020
%==========================================================================

function H = whisart_entropy( parameters)
m = parameters.q;
L = parameters.L;
g = 0;
for ii=1:parameters.q
    aux = gamma((parameters.L-(ii-1))/2);
    if(isinf(aux))
        aux = 0;
    end
    g = g+aux;
end
H = ( m*(m-1)/2 )*log(pi) + (m*m*log(L)) + (m*log(det(parameters.cov))) + g;
end

