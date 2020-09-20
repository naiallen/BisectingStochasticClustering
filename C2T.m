%%=========================================================================
% Brief: Convert from covariance matrix to coherence matrix
%   Input:
%       cov: covariance matrix 3x3
%   Output
%       coh: coherence matrix 3x3
% Autor: Naiallen Carvalho
% Computação Aplicada
% Instituto Nacional de Pesquisas Espaciais - INPE
% 2020
%==========================================================================
function [ coh ] = C2T( cov )
C11 = cov(:,1);
C12_re = real(cov(:,2));
C12_im = imag(cov(:,2));
C13_re = real(cov(:,3));
C13_im = imag(cov(:,3));
C22 =  cov(:,5);
C23_re = real(cov(:,6));
C23_im = imag(cov(:,6));
C33 = cov(:,9);

coh = zeros(size(cov));

coh(:,1) = (C11 + 2 * C13_re + C33) / 2;
coh(:,2) = (C11 - C33) / 2 - 1i*C13_im;
coh(:,3) = (C12_re + C23_re) / sqrt(2) + 1i*(C12_im - C23_im) / sqrt(2);
coh(:,4) = conj(coh(:,2));
coh(:,5) = (C11 - 2 * C13_re + C33) / 2;
coh(:,6) = (C12_re - C23_re) / sqrt(2) + 1i*(C12_im + C23_im) / sqrt(2);
coh(:,7) = conj(coh(:,3));
coh(:,8) = conj(coh(:,6));
coh(:,9) = C22;
end

