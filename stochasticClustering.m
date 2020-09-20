%%=========================================================================
% Brief: Run Stochastic Clustering
%   Input:
%       input_im
%       distance
%       n_classes
%       seeds
%       iteration
%       L
%       Root_ent
%   Output
%       index
% Autor: Naiallen Carvalho
% Computação Aplicada
% Instituto Nacional de Pesquisas Espaciais - INPE
% 2020
%==========================================================================
function [index, output_im] = stochasticClustering(input_im, pdf_name, distance, n_classes, seeds, iteration, L)

ampl_max = [max(input_im(:,1)) max(input_im(:,5)) max(input_im(:,9))];
ampl_min = [min(input_im(:,1)) min(input_im(:,5)) min(input_im(:,9))];
ampl_max(isnan(ampl_max)) = 0;
ampl_min(isnan(ampl_min)) = 0;
error_th = mean(abs(ampl_max - ampl_min))*0.05;
[n_row, n_col, n_bands] = size(input_im);
index =1;
counter = 1;
data_total = im2vec( input_im);
% covariance = reshape(intrisic_mean( data_total ),3 ,3);
% par.L = L;
% par.q =  3;
% par.cov = covariance;
% H = abs(whisart_entropy( par, n_bands) / Root_ent);
% w_size = (size(data_total,1) /(n_row*n_col));
step = floor(L/2);
while (index <= iteration)
    output_im = ones(n_row, n_col, n_classes)*nan;
    for ii = 1:n_row          %linha
        for jj = 1:n_col      %coluna
            if isnan(input_im(ii, jj,:))
                continue;
            end
            
            %% Get covariance
            dist = zeros(2,1);
            if (strcmp(pdf_name, 'wishart'))
                if (n_bands == 9)
                    cov1 = reshape(input_im(ii, jj,:), 3,3);
                elseif(n_bands == 18)
                    cov1 = zeros(6, 6);
                    cov1(1:3, 1:3) = reshape(input_im(ii, jj, 1:9), 3,3);
                    cov1(4:6, 4:6) = reshape(input_im(ii, jj, 10:18), 3,3);
                end
                for kk = 1:n_classes
                    cov2 = seeds{kk};
                    dist(kk, 1) = abs( wishart_stochastic_distance(distance, cov1, cov2, L) );
                end
            elseif(strcmp(pdf_name, 'gaussian'))
                if (ii < step+1)||(jj < step+1)
                    continue;
                end
                if (ii > (n_row - step)) || (jj > (n_col - step))
                    continue;
                end
                temp_data = im2vec(input_im(ii-step:ii+step, jj-step:jj+step, :));
                cov1 = get_covariance(temp_data);
                m1 = mean(temp_data);
                for kk = 1:n_classes
                    cov2 = seeds{kk}.cov;
                    m2 = seeds{kk}.m;
                    dist(kk, 1) = abs( gaussian_stochastic_distance(distance, m1, m2, cov1, cov2) );
                end
            end
            [~, pos] = min(dist);
            output_im(ii, jj, pos) = 1;
            
        end
    end
    
    %% Geet new seeds
    ampl_prev = [];
    ampl_actual = [];
    
    for kk=1:n_classes
        im = zeros(n_row, n_col, n_bands);
        for nb = 1:n_bands
            im(:, :, nb) = output_im(:,:, kk).*(input_im(:,:,nb));
        end
        clear data
        data = im2vec( im );
        if (strcmp(pdf_name, 'wishart'))            
            if (n_bands == 9)
                ampl_prev = [ampl_prev seeds{kk}(1,1) seeds{kk}(2,2) seeds{kk}(3,3)];
            elseif(n_bands == 18)
                ampl_prev = [ampl_prev seeds{kk}(1,1) seeds{kk}(2,2) seeds{kk}(3,3) seeds{kk}(4,4) seeds{kk}(5,5) seeds{kk}(6,6)];
            end

            if((~isempty(data)) && (size(data,1) > 5))
                seeds{kk} = intrisic_mean( data );
            else
                s = size(data_total, 1);
                index1 = randi(s);
                m_cov = data_total(index1, :);
                if (n_bands == 9)
                    seeds{kk} = reshape(m_cov,3,3)';
                elseif(n_bands == 18)
                    cov = zeros(6, 6);
                    cov(1:3, 1:3) = reshape(m_cov(1:9),3,3);
                    cov(4:6, 4:6) = reshape(m_cov(10:18),3,3);
                    seeds{kk} = cov;
                end
                counter = counter+1;
            end
            if (n_bands == 9)
                ampl_actual = [ampl_prev seeds{kk}(1,1) seeds{kk}(2,2) seeds{kk}(3,3)];
            elseif(n_bands == 18)
                ampl_actual = [ampl_prev seeds{kk}(1,1) seeds{kk}(2,2) seeds{kk}(3,3) seeds{kk}(4,4) seeds{kk}(5,5) seeds{kk}(6,6)];
            end
        elseif(strcmp(pdf_name, 'gaussian'))
            ampl_prev = [ampl_actual seeds{kk}.m(1) seeds{kk}.m(2) seeds{kk}.m(3)];
            seeds{kk}.cov = get_covariance( data );
            seeds{kk}.m = mean( data );
            ampl_actual = [ampl_actual seeds{kk}.m(1) seeds{kk}.m(2) seeds{kk}.m(3)];
        end
    end
    %% Compute RSME
    error = sqrt( mean((abs(ampl_prev.^2 - ampl_actual.^2))));
    if((error < error_th)||...
            (index >= 20)||...
            (counter >= 3) )
        break;
    end
    %     error_prev = error;
    index=index+1;
end
index=index-1;



end

