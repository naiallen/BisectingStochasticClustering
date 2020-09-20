%%=========================================================================
% Brief: Compute stochastic distance
%   Input:
%       distance_name
%       cov1
%       cov2
%       L
%   Output
%       value
% Autor: Naiallen Carvalho
% Computação Aplicada
% Instituto Nacional de Pesquisas Espaciais - INPE
% 2020
%==========================================================================
%
function value = gaussian_stochastic_distance(distance_name, m1, m2, cov1, cov2)
q = size(cov1,1);
cov1 = cov1+eye(q)*0.001;
cov2 = cov2+eye(q)*0.001;
media = m1.' - m2.';
cov = cov1 + cov2;
%% Bhattacharyya
if ( strcmp(distance_name, 'Bhattacharyya') )
    value = (1/8)*( media.'*inv(cov*0.5)*media ) + 0.5*log(det(0.5*cov)/sqrt(det(cov1) * det(cov2)) );
    %% Kullback_Leibler
elseif(strcmp(distance_name, 'Kullback-Leibler') )
    value = 0.5*( media.'*(inv(cov1) + inv(cov2))*media ) + 0.5*(trace(inv(cov1)*cov2+inv(cov2)*cov1)-(2*size(cov2,1)));
    %% Euclidean Distance
elseif strcmp(distance_name, 'Euclidean')
    value = sqrt(sum((mu_1 - mu_2 ) .^ 2));
end
end