%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DPC_tik solves the least squared solutions for phase and absorption     %
% with Tikhonov regularization                                            %
%              ||        ||2        || ||2      || ||2                    %
% Problem: min ||fIDPC-Ax||  + alpha||a|| + beta||p||                     %
%           x  ||        ||         || ||       || ||                     %
%                                                                         %
%            H    |M11 M12|       |a|                                     %
%           A  A= |M21 M22| , x = |p|                                     %
% a: F{amplitude}, p: F{phase}                                            %
%                                                                         %
% Inputs:                                                                 %
%        fIDPC        : Fourier spectrum of normalized measurements       %
%        source       : illumination patterns                             %
%        pupil        : pupil function defined by the NA of imaging system%
%        reg_tik      : regularization parameters (alpha, beta)           %
% Outputs:                                                                %
%        absorption   : absorption of the complex field                   %
%        phase        : phase of the complex field                        %
%                                                                         %
% Copyright (C) 2018 Michael Chen                                         %
%                                                                         %
% This program is free software: you can redistribute it and/or modify    %
% it under the terms of the GNU General Public License as published by    %
% the Free Software Foundation, either version 3 of the License, or       %
% (at your option) any later version.                                     %
%                                                                         %
% This program is distributed in the hope that it will be useful,         %
% but WITHOUT ANY WARRANTY; without even the implied warranty of          %
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           %
% GNU General Public License for more details.                            %
%                                                                         % 
% You should have received a copy of the GNU General Public License       %
% along with this program.  If not, see <http://www.gnu.org/licenses/>.   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [absorption, phase] = DPC_tik(fIDPC, source, pupil, reg_tik)
    IF              = @(x) ifft2(x);
    Hi              = zeros(size(source));
    Hr              = zeros(size(source));    
    
    % evaluate analytical transfer functions
    for source_index = 1:size(source, 3)
        [Hi(:, :, source_index),...
         Hr(:, :, source_index)] = genTransferFunction(source(:, :, source_index), pupil);
    end
    
    % matrix pseudo inverse
    M11             = sum(abs(Hr).^2, 3) + reg_tik(1);
    M12             = sum(conj(Hr).*Hi, 3);
    M21             = sum(conj(Hi).*Hr, 3);
    M22             = sum(abs(Hi).^2, 3) + reg_tik(2);
    denominator     = M11.*M22 - M12.*M21;
    I1              = sum(fIDPC.*conj(Hr), 3);
    I2              = sum(fIDPC.*conj(Hi), 3);
    absorption      = real(IF((I1.*M22-I2.*M12)./denominator));
    phase           = real(IF((I2.*M11-I1.*M21)./denominator));   
end

