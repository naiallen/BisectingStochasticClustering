%%=========================================================================
% Brief: Plot PolSAR image
%   Input: 
%       PolSAR image
%   Output
%       plot
% Autor: Naiallen Carvalho
% Computação Aplicada
% Instituto Nacional de Pesquisas Espaciais - INPE
% 2020
%==========================================================================
function plotPolSARimage( input_im )
%% Get size
[~, ~, n_bands]=size(input_im);

%% Get Image
if (n_bands == 1)
    HH = 10*log10(abs(input_im(:,:, 1)));
    HV = 10*log10(abs(input_im(:,:, 1)));
    VV = 10*log10(abs(input_im(:,:, 1)));
elseif (n_bands == 3)
    HH = 10*log10(abs(input_im(:,:, 1)));
    HV = 10*log10(abs(input_im(:,:, 2)));
    VV = 10*log10(abs(input_im(:,:, 3)));
elseif(n_bands == 9) || (n_bands == 18)
    HH = 10*log10(abs(input_im(:,:, 1)));
    HV = 10*log10(abs(input_im(:,:, 5)));
    VV = 10*log10(abs(input_im(:,:, 9)));
end

%%
color = zeros(size(HH, 1), size(HH, 2), 3);
XMIN = min([min(HH(:)) min(HV(:)) min(VV(:))]);
XMAX = max([max(HH(:)) max(HV(:)) max(VV(:))]);
color(:,:,1) = (HH-XMIN)/(XMAX-XMIN);
color(:,:,2) = (HV-XMIN)/(XMAX-XMIN);
color(:,:,3) = (VV-XMIN)/(XMAX-XMIN);
imagesc(color); %Mostra a imagem em escala
% axis equal
camroll(90)
% set(gca, 'xdir', 'reverse');
axis off
end

