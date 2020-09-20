function [index, seed, data1, data2] = get_seed_EM_gaussian(imacov_vec, q)

[index, data1, data2] = EM_gaussian( imacov_vec, q );

if (isempty(data1))
    %% Get Seed 2----------------------------------------------------------
    seed{2}.m = mean(data2);
    seed{2}.cov = get_covariance(data2);
    
    %% Get Seed 1----------------------------------------------------------
    seed{1} = seed{2};
    data1 = data2(1:3,:);
elseif (isempty(data2))
    %% Get Seed 1----------------------------------------------------------
    seed{1}.m = mean(data1);
    seed{1}.cov = get_covariance(data1);
    
    %% Get Seed 2----------------------------------------------------------
    seed{2}= seed{1};
    data2 = data1(1:3,:);
else
    %% Get Seed 1----------------------------------------------------------
    seed{1}.m = mean(data1);
    seed{1}.cov = get_covariance(data1);
    
    %% Get Seed 2----------------------------------------------------------
    seed{2}.m = mean(data2);
    seed{2}.cov = get_covariance(data2);
end

end

