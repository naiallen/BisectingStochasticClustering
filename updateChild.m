function tree = updateChild(tree, next_level, child_index, parent_index, data_cov, L, q,  n_row, n_col, initial_parameter_algo, image_type)

if (size(data_cov,1) <= 5)
    tree{next_level, child_index}.EntropyGain = 0;
    return
end

%% Get seeds
if (strcmp(image_type, 'PolSAR'))
    if(strcmp(initial_parameter_algo, 'RPDDP'))
        [ch_seeds, data1, data2] = get_seed_RPDDP(data_cov);
    elseif (strcmp(initial_parameter_algo, 'EM'))
        [~, ch_seeds] = get_seed_EM_wishart(data_cov,  L, q);
    end
elseif(strcmp(image_type, 'Optical'))
    if(strcmp(initial_parameter_algo, 'RPDDP'))
        [ch_seeds, data1, data2] = get_seeds_pca_real(data_cov);
    elseif (strcmp(initial_parameter_algo, 'EM'))
        [~, ch_seeds, data1, data2] = get_seed_EM_gaussian(data_cov, q);
    end
    
end

tree{next_level, child_index}.Seeds = ch_seeds;
% tree{next_level, child_index}.SeedsData = {data1 data2};

%% Get parent parameter
if (strcmp(image_type, 'PolSAR'))
    m_cov = intrisic_mean( data_cov );
    covariance = m_cov;% average_covariance(data);
    par.L = L;
    par.q =  q;
    par.cov = covariance;
    tree{next_level, child_index}.Parameter = par;
    tree{next_level, child_index}.Entropy = whisart_entropy(par);
elseif(strcmp(image_type, 'Optical'))
    par.q =  q;
    par.m =  mean(data_cov);
    par.cov = get_covariance(data_cov);
    entropy1 = gaussian_entropy( par);
    tree{next_level, child_index}.Entropy = gaussian_entropy(par);
end

%% Get entropy


%% Get information gain parent
tree{next_level, child_index}.EntropyGainParent = abs( tree{next_level-1, parent_index}.Entropy - tree{next_level, child_index}.Entropy );
if(isnan(tree{next_level, child_index}.EntropyGainParent))
    test = 0;
end
if (tree{next_level, child_index}.EntropyGainParent < 1e-1)
    tree{next_level, child_index}.EntropyGainParent = 0;
end


%% Compute entropy gain - child 1
if (strcmp(image_type, 'PolSAR'))
    par.L = L;
    par.q =  q;
    par.cov = ch_seeds{1};
    entropy1 = whisart_entropy(par);
elseif(strcmp(image_type, 'Optical'))
    par.q =  q;
    par.m = ch_seeds{1}.m;
    par.cov = ch_seeds{1}.cov;
    entropy1 = gaussian_entropy( par);
end        


%% Get size weigth - child 1
sw1 = size(data1, 1)/(n_row*n_col);

%% Get entropy gain - child 1
HG1 = abs((tree{next_level, child_index}.Entropy-entropy1)*sw1);

%% Compute entropy gain - child 2
if (strcmp(image_type, 'PolSAR'))
    par.L = L;
    par.q =  q;
    par.cov = ch_seeds{2};
    entropy2 = whisart_entropy(par);
elseif(strcmp(image_type, 'Optical'))
    par.q =  q;
    par.m = ch_seeds{2}.m;
    par.cov = ch_seeds{2}.cov;
    entropy2 = gaussian_entropy( par);
end      


%% Get size weigth - child 2
sw2 = size(data2, 1)/(n_row*n_col);

%% Get entropy gain - child 2
HG2 = abs((tree{next_level, child_index}.Entropy-entropy2)*sw2);

HG = mean([HG1 HG2]);
if ( (HG < 0.01))%  && ( (size(data_cov,1)/(n_row*n_col) ) < 0.15 ) )
    HG = 0;
    % % % elseif (HG < 0.6)
    % % %     HG = 0;
end
tree{next_level, child_index}.EntropyGain = HG;
end

