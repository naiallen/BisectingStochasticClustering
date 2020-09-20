%% ========================================================================
% Brief: Bisecting Stochasti Clustering
%   Input:
%       Covariance Image
%   Output
%       Dendrogram
%       Classified Image
% Autor: Naiallen Carvalho
% Computação Aplicada
% Instituto Nacional de Pesquisas Espaciais - INPE
% 2020
%==========================================================================

clear
close all
clc

%% 0. defines==============================================================
TRUE = 1;
FALSE = 0;
%==========================================================================

%% 1. Initial configuration================================================
level_max = 30;
use_size_weigth = TRUE;
enable_plan_H_alpha = FALSE;
initial_parameter_algo = 'EM'; %'RPDDP'; %
image_type =  'Optical'; %'PolSAR'
L = 5; %Number of Looks;
distance_name = 'Bhattacharyya'; % 'Hellinger';
folder = 'C:\Users\ncarval1\Documents\PersonalCodes\Matlab\Doutorado\Imagens\CovarianceImages\';
folderTosave = 'C:\Users\ncarval1\Documents\PersonalCodes\Matlab\Doutorado\Imagens\CovarianceImages\test\';

%==========================================================================

%% 3. Load image===========================================================
if (strcmp(image_type, 'PolSAR'))
    im_filename = 'SIRC_Lband_L5.dat'; %'ALOS_SanFrancisco_L5.dat'; %
    hdr_filename = 'SIRC_Lband_L5.hdr'; %'ALOS_SanFrancisco_L5.hdr'; %
    imcov = openPolSARimage(folder, im_filename, hdr_filename);
    [n_row, n_col, n_band] = size(imcov);
    figure;
    plotPolSARimage( imcov );
    if (n_band == 9)
        q = 3;
    elseif(n_band == 18)
        q = 6;
    end
    
    pdf_name = 'wishart';
elseif (strcmp(image_type, 'Optical'))
    im_filename = 'san_francisco.jpg';
    address = fullfile(folder, im_filename);
    imcov = double(imread(address))/255;
    figure;
    imagesc(imcov);
    axis equal
    [n_row, n_col, n_band] = size(imcov);
    q = 3;
    pdf_name = 'gaussian';
end
%==========================================================================

%% 4. Enable Plan H-alpha
if (n_band == 9) && (strcmp(image_type, 'PolSAR'))
    enable_plan_H_alpha = TRUE;
end
%==========================================================================


%% 5. Compute coherence matrix image and compute entropy and alpha angle===
if(enable_plan_H_alpha)
    % 3.1 Get coherence matrix
    imcov_vec = im2vec(imcov);
    imcoh_vec = C2T( imcov_vec );
    imcoh = reshape(imcoh_vec, size(imcov));
    figure;
    plotPolSARimage( imcoh );
    
    % 3.1 Compute the entropy and alpha angle
    [entropia, alpha] = planHalpha(imcoh_vec);
    figure;
    plotPlanHalpha( entropia, alpha )
end
%==========================================================================

%% 6.  Algoritimo hierárquico==============================================
tree = [];
parent_id = 1;
im_dummy = zeros(n_row, n_col);
numinte = zeros(level_max,1);
for level = 1:level_max
    next_level = level+1;
    %% Split levels
    if level ==  1
        %% Populate tree---------------------------------------------------
        tree = populate_tree(tree, parent_id, level);
        parent_index = 1;
        
        %% Vectorize image-------------------------------------------------
        imcov_vec = im2vec(imcov);
        
        %% Get parent parameter--------------------------------------------
        % PolSAR image
        if (strcmp(image_type, 'PolSAR'))
            par.L = L;
            par.q =  q;
            par.cov = intrisic_mean( imcov_vec );
        elseif(strcmp(image_type, 'Optical'))
            par.q =  q;
            par.m =  mean(imcov_vec);
            par.cov = get_covariance( imcov_vec );
        end
        
        tree{parent_index}.Parameter = par;
        
        %% Get seeds
        if (strcmp(image_type, 'PolSAR'))
            if(strcmp(initial_parameter_algo, 'RPDDP'))
                [seeds, ~, ~] = get_seed_RPDDP(imcov_vec);
            elseif (strcmp(initial_parameter_algo, 'EM'))
                [indexEM, seeds] = get_seed_EM_wishart(imcov_vec,  L, q);
            end
        elseif(strcmp(image_type, 'Optical'))
            if(strcmp(initial_parameter_algo, 'RPDDP'))
                [seeds, ~, ~] = get_seeds_pca_real(imcov_vec);    
            elseif (strcmp(initial_parameter_algo, 'EM'))
                [indexEM, seeds, ~, ~] = get_seed_EM_gaussian(imcov_vec, q);
            end
            
                    
        end
        tree{parent_index}.Seeds = seeds;
        
        %% Get entropy parent
        if (strcmp(image_type, 'PolSAR'))
            tree{parent_index}.Entropy = whisart_entropy( par);
        elseif(strcmp(image_type, 'Optical'))
            tree{parent_index}.Entropy = gaussian_entropy( par);
        end
        
        %% Populate dummy image
        im_dummy = ones(n_row, n_col);
        
        %% Get child ID
        child1_id = tree{parent_index}.ChildsID(1);
        child2_id = tree{parent_index}.ChildsID(2);
        
        child1_index = getTreeIndex( tree, next_level, child1_id);
        child2_index = getTreeIndex( tree, next_level, child2_id);
    else
        
        %% Choose child to open
        n_leaf = size(tree, 2);
        HG = zeros(n_leaf, 1);
        
        for ii = 1:n_leaf
            HG(ii) = tree{level, ii}.EntropyGain;
        end
        [~, pos] = max(HG);
        parent_index = pos;
        
        %% Populate tree
        parent_id = tree{level, parent_index}.ID;
        tree = setChildID(tree, parent_id, level);
        tree = populate_tree(tree, parent_id, level);
        
        child1_id = tree{level, parent_index}.ChildsID(1);
        child2_id = tree{level, parent_index}.ChildsID(2);
        
        child1_index = getTreeIndex( tree, next_level, child1_id);
        child2_index = getTreeIndex( tree, next_level, child2_id);
        %% Get seeds
        seeds = tree{level, parent_index}.Seeds;
    end
    
    
    im_aux = im_dummy;
    im_aux(im_aux ~= parent_id) = NaN;
    im_test = ones(n_row, n_col, n_band);
    for nb = 1:n_band
        im_test(:,:,nb) = im_aux.*(imcov(:,:,nb));
    end
    if(level == 14)
        test = 0;
    end
    %% Plot parent image---------------------------------------------------
    h1 = figure('Position', [10 10 1400 700]);
    subplot(242)
    if (strcmp(image_type, 'PolSAR'))
        plotPolSARimage(im_test);
    elseif(strcmp(image_type, 'Optical'))
        imagesc(im_test);
    end
    t1 = '';
    title({strcat(strcat(strcat('Parent - Level: ', num2str(level)), ' - ID:'), num2str(tree{level, parent_index}.ID)), t1}, 'color', tree{level, parent_index}.Color);
    
    
    %% Run Stochastic Clustering-------------------------------------------
    tic
    [numinte(level), output_im] = stochasticClustering(im_test, pdf_name, distance_name, 2, seeds, 20, L);
    %     numinte(level) = numinte(level)+indexEM;
    toc
    %----------------------------------------------------------------------
    
    %% Get Image ----------------------------------------------------------
    im1_cov = zeros(n_row, n_col, n_band);
    d = output_im(:,:,1);
    d1 = d(:);
    d = output_im(:,:,2);
    d2 = d(:);
    for nb = 1:n_band
        im1_cov(:,:,nb) = output_im(:,:,1).*(imcov(:,:,nb));
    end
    
    im2_cov = zeros(n_row, n_col, n_band);
    for nb = 1:n_band
        im2_cov(:,:,nb) = output_im(:,:,2).*(imcov(:,:,nb));
    end
    %----------------------------------------------------------------------
    
    %% Get vector data ----------------------------------------------------
    data1_cov = im2vec( im1_cov);
    data2_cov = im2vec( im2_cov);
    %----------------------------------------------------------------------
    
    subplot(246)
    if (strcmp(image_type, 'PolSAR'))
        plotPolSARimage(im1_cov);
    elseif(strcmp(image_type, 'Optical'))
        imagesc(im1_cov);
    end
    
    
    subplot(247)
    if (strcmp(image_type, 'PolSAR'))
        plotPolSARimage( im2_cov );
    elseif(strcmp(image_type, 'Optical'))
        imagesc(im2_cov);
    end   
    
    if(enable_plan_H_alpha)
        hold on
        data1 = [entropia(~isnan(d1)) alpha(~isnan(d1))];
        data2 = [entropia(~isnan(d2)) alpha(~isnan(d2))];
        data = [data1; data2];
        
        subplot(243)
        plotPlanHalpha( data(:,1), data(:,2));
        
        ent1 = 0;
        ent2 = 0;
        alp1 = 0;
        alp2 = 0;
        if (isempty(data1) || isempty(data2))
            error_entropy = 0;
            error_alpha = 0;
        else
            ent1 = mean(data1(:,1));
            alp1 = mean(data1(:,2));
            
            ent2 = mean(data2(:,1));
            alp2 = mean(data2(:,2));
            
            error_entropy = round(abs(ent1 - ent2),2);
            error_alpha = round(abs(alp1 - alp2), 2);
        end
        
        subplot(244)
        plotHalphaOnly;
        hold on
        plot(  ent1,  alp1, 'o', 'LineWidth',2, 'color', tree{next_level, child1_index}.Color );
        hold on
        plot(  ent2,  alp2, 'o', 'LineWidth',2, 'color', tree{next_level, child2_index}.Color );
        plot(  mean(data(:,1)),  mean(data(:,2)), 'o', 'LineWidth',2, 'color', tree{level, parent_index}.Color );
        title({strcat('Entropy Error: ', num2str(error_entropy)),strcat('Alpha Angle Error: ', num2str(error_alpha))})
        
        subplot(245)
        plotPlanHalpha( data1(:,1), data1(:,2) )
        
        subplot(248)
        plotPlanHalpha(data2(:,1), data2(:,2))
    end
    
    
    %% update image dummy
    im_dummy(output_im(:,:,1)== 1) = child1_id;
    im_dummy(output_im(:,:,2)== 1) = child2_id;
    
    %% Get childs indexes
    child_id = tree{level, parent_index}.ChildsID(1);
    child1_index = getTreeIndex( tree, level+1, child_id);
    child_id = tree{level, parent_index}.ChildsID(2);
    child2_index = getTreeIndex( tree, level+1, child_id);
    
    %% Compute the child 1 seeds and parameters----------------------------
    w1_d1 = size(data1_cov,1)/(n_row*n_col);
    w1_d2 = size(data2_cov,1)/(n_row*n_col);
    w2_d1 = size(data1_cov,1)/(size(data1_cov,1)+size(data2_cov,1));
    w2_d2 = size(data2_cov,1)/(size(data1_cov,1)+size(data2_cov,1));
    
    next_level = level+1;
    
    
    if (isempty(data1_cov) || isempty(data2_cov) ||...
            (w2_d1 < 0.01) || (w2_d2 < 0.01))
        %% If there is no division, the branch stops
        tree{next_level, child1_index}.EntropyGain = 0;
        tree{next_level, child2_index}.EntropyGain = 0;
    else
        %% Update child node
        if (w1_d1 < 0.005)
            %% If there is no division, the branch stops
            tree{next_level, child1_index}.EntropyGain = 0;
        else
            tree = updateChild(tree, next_level, child1_index, parent_index, data1_cov, L, q, n_row, n_col, initial_parameter_algo, image_type);
        end
        if (w1_d2 < 0.005)
            tree{next_level, child2_index}.EntropyGain = 0;
        else
            tree = updateChild(tree, next_level, child2_index, parent_index, data2_cov, L, q,  n_row, n_col, initial_parameter_algo, image_type);
        end
    end
    
    subplot(246)
    t1 = strcat( 'HG (Parent - Child):', num2str(round(tree{next_level, child1_index}.EntropyGainParent, 2)) );
    title({strcat(strcat(strcat('Child 1 - Level: ', num2str(next_level)), ' - ID:'), num2str(tree{next_level, child1_index}.ID)), t1}, 'color', tree{next_level, child1_index}.Color);
    
    subplot(247)
    t1 = strcat( 'HG (Parent - Child):', num2str(round(tree{next_level, child2_index}.EntropyGainParent, 2)) );
    title({strcat(strcat(strcat('Child 2 - Level: ', num2str(next_level)), ' - ID:'), num2str(tree{next_level, child2_index}.ID)), t1}, 'color', tree{next_level, child2_index}.Color);
    
    
    saveas(h1, strcat(folderTosave, strcat(strcat('level_', num2str(level)), '.png')));
    HG = 0;
    for index =1:(size(tree, 2))
        if (~isempty(tree{next_level, index}))
            HG = HG + tree{next_level, index}.EntropyGain;
        end
    end
    if (HG <= 0)
        break;
    end
end
h4  = figure('Position', [10 10 2000 700]);
plotTree(tree);
saveas(h4, strcat(folderTosave, 'tree.png'));

h3 = figure;
plot(numinte)
xlabel('Level')
ylabel('Number of iteration');
grid on
saveas(h3, strcat(folderTosave, 'numberofiteration.png'));

h2  = figure();
imagesc (im_dummy);
colormap gray
camroll(-90);
set(gca, 'xdir', 'reverse');
axis equal
imwrite(uint8(im_dummy), strcat(folderTosave, 'im_classifgray.png'))


