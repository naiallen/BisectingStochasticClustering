%%=========================================================================
% Brief: Get tree index
%   Input:
%       tree
%       level
%       ID
%   Output
%       index
% Autor: Naiallen Carvalho
% Computação Aplicada
% Instituto Nacional de Pesquisas Espaciais - INPE
% 2020
%==========================================================================
function index = getTreeIndex( tree,level, ID)
%% look for ID
for ii = 1:(size(tree, 2))
    if (~isempty(tree{level, ii}))
        if (tree{level, ii}.ID == ID)
            index = ii;
            break;
        end
    end
end
end

