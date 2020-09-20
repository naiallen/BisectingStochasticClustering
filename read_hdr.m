%%=========================================================================
% Brief: Read PolSAR image header
%   Input: 
%       hdr_filename: Image SAR data header address
%   Output
%       n_row: numbe of row
%       n_col: number of colunm
%       n_band: number of bands
%       format: image format (float or double)
%       interleave: bsq (Band-Sequential) or bip (Band-Interleaved-by-Pixel)
%       byte_order: 0 ( native format)  and 1 (IEEE format)
% Autor: Naiallen Carvalho
% Computação Aplicada
% Instituto Nacional de Pesquisas Espaciais - INPE
% 2020
%==========================================================================
function [ n_row, n_col, n_band, format, interleave, byte_order] = read_hdr( file_addr )
fileID = fopen(file_addr);
tline = fgetl(fileID);
while ischar(tline)
    str = strcat(tline);
    if(contains( str , 'samples' ))
        n_row = str2double(regexprep(str, 'samples =',''));
    end
    
    if(contains( str , 'lines' ))
        n_col = str2double(regexprep(str, 'lines   =',''));
    end
    
    if(contains( str , 'bands' ))
        n_band = str2double(regexprep(str, 'bands   =',''));
    end
    
    if(contains( str , 'data type' ))
        type = str2double(regexprep(str, 'data type =',''));
        if(type == 6)
            format = 'float';
        elseif(type == 9)
            format = 'double';
        end
    end
    
    if(contains( str , 'interleave ' ))
        interleave  = (regexprep(str, 'interleave =',''));
    end
    
    if(contains(str , 'byte order ' ))
        byte_order  = str2double(regexprep(str, 'byte order = ',''));
    end

    tline = fgetl(fileID);
end
fclose(fileID);
end