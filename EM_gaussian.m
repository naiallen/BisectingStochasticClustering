function [index, data1, data2] = EM_gaussian( imacov_vec, q )
s = size(imacov_vec, 1);

index1 = randi(s);
index2 = randi(s);

%% Get parameter 1---------------------------------------------------------
p{1}.weight = 0.5;
p{1}.q = q;
p{1}.m = imacov_vec(index1, :);
p{1}.cov =get_covariance(imacov_vec(index1-2:index1+2,:));

%% Get parameter 2---------------------------------------------------------
p{2}.weight = 0.5;
p{2}.q = q;
p{2}.m = imacov_vec(index2, :);
p{2}.cov = get_covariance(imacov_vec(index2-2:index2+2,:));

%%
epsilon = 1e-2; % precision

error = 2;
index = 0;
step = ceil(size(imacov_vec,1)*0.05);
while (error > epsilon)
    %% Expectation---------------------------------------------------------
    data1 = [];
    data2 = [];
    for ii = 1:step:size(imacov_vec,1)
        
        x = imacov_vec(ii, :);
        p_cluster1 = log(p{1}.weight) - getLogGaussian( x, p{1} );
        p_cluster2 = log(p{2}.weight) - getLogGaussian( x, p{2} );
                
        if abs(p_cluster1) > abs(p_cluster2)
            data1 = [data1; imacov_vec(ii, :)];
        else
            data2 = [data2; imacov_vec(ii, :)];
        end
        
    end
    
    if ( isempty(data1) || isempty(data2) )
        break
    end
    %----------------------------------------------------------------------
    
    %% Maximization--------------------------------------------------------
    %weigths
    p{1}.weight = 0.5;%size(data1,1) / (size(imacov_vec,1)/10);
    p{2}.weight = 1 - p{1}.weight;
    
    prev_temp = [p{1}.m p{2}.m]';
    if((~isempty(data1)) && (size(data1,1)>1))
        p{1}.m = mean( data1 );
        p{1}.cov = get_covariance( data1 );
    end
    
    if((~isempty(data2)) && (size(data2,1)>1))
        p{2}.m = mean( data2 );
        p{2}.cov = get_covariance( data2 );
    end
    
    actual_temp = [p{1}.m p{2}.m]';
    temp_num = sum( ((prev_temp - actual_temp) .^ 2)' );
    tem_den = sum( (prev_temp.^2)' );
    error_v = temp_num./tem_den;
    error = max(error_v);
    index = index+1;
    if(index > 10)
        break;
    end
end

end

