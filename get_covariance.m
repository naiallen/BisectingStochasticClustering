%% ========================================================================
% Brief: Compute covariance
%   Input:
%       parameters
%   Output
%       H: entropy
% Autor: Naiallen Carvalho
% Computação Aplicada
% Instituto Nacional de Pesquisas Espaciais - INPE
% 2020
%==========================================================================

function sigma = get_covariance(data_vec)
    n_bands = size(data_vec, 2);
    if (n_bands > 1 && size(data_vec, 1) > 1)
        m =  mean( data_vec );
        for ii = 1:n_bands
            data_vec(:, ii) = data_vec(:, ii) - m(ii);
        end
        sigma = ( 1/(size(data_vec, 1)-1) )*(data_vec'*data_vec);
    else
        sigma = nan;
    end
end

