%%=========================================================================
% Brief: Compute intrisic mean
%   Input: 
%       imacov_vec
%   Output
%       cov
% Autor: Naiallen Carvalho
% Computação Aplicada
% Instituto Nacional de Pesquisas Espaciais - INPE
% 2020
%==========================================================================

function cov = intrisic_mean( imacov_vec )
%% Intrisic mean
n_bands = size(imacov_vec, 2);

N_iteration = 0;
if (n_bands == 9)
    N_iteration = 1;
    cov = zeros(3, 3);
elseif(n_bands==18)
    N_iteration = 2;
    cov = zeros(6, 6);
end

for jj = 1:N_iteration
    
    if (n_bands == 9)
        imacov_vec_temp = imacov_vec;
    elseif(n_bands==18)
        if(jj == 1)
            imacov_vec_temp = imacov_vec(:, 1:9);
        elseif(jj == 2)
            imacov_vec_temp = imacov_vec(:, 10:18);
        end
    end
    
    N = size(imacov_vec_temp, 1);
    m_cov_temp = sum(imacov_vec_temp)/N;
    m_cov = reshape(m_cov_temp,3,3);
    error = 100;
    index = 0;
    while(error > 0.01)
        n_bands_temp = 9;
        delta_m = zeros(sqrt(n_bands_temp), sqrt(n_bands_temp));
        for ii =1:N
            A = reshape(imacov_vec_temp(ii,:),3,3);
            if ((sum(eig(A) < 0)) > 0)
                continue
            end
            delta_m = delta_m + (logm(A));
        end
        delta_m = delta_m/N;
        m_cov_prev = m_cov;
        m_cov = expm(delta_m);
        error = abs( trace(m_cov_prev) - trace(m_cov) );
        index = index+1;
        if (index > 5)
            break;
        end
    end
    if(trace(m_cov) == 3)
        m_cov = average_covariance(imacov_vec_temp);
    end
    
    if (n_bands == 9)
        cov = m_cov;
    elseif(n_bands==18)
        if(jj == 1)
            cov(1:3, 1:3) = m_cov;
        elseif(jj == 2)
            cov(4:6, 4:6) = m_cov;
        end
    end    
end
