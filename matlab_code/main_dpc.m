%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main_dpc -                                                              %
% Main file for DPC absorption and phase recovery                         %
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
clear; close all;
tic;
set(0, 'DefaultFigureWindowStyle', 'docked');
addpath('dpc_functions');
F             = @(x) fft2(x);
IF            = @(x) ifft2(x);

%% load data
% load DPC images, subtract DC, and normalized by the total energy

data_path     = ['..', filesep, 'sample_data']; %PUT YOUR DATA PATH HERE

image_list    = dir([data_path, filesep, '*.tif']);
for image_index = 1:numel(image_list)
   image_load = double(imread([data_path, filesep, image_list(image_index).name]));
   if image_index==1
       IDPC                    = zeros(size(image_load, 1), size(image_load, 2), numel(image_list));
   end
   IDPC(:, :, image_index) = image_load/mean2(image_load)-1;
end

%% system parameters
dim           = [size(IDPC, 1), size(IDPC, 2)]; % image size
sigma         = 1.0;                            % partial coherence factor
na            = 0.40;                           % numerical aperture of the imaging system
na_illum      = sigma*na;                       % numerical aperture of the illumination
magnification = 20*2;                           % magnification of the imaging system
lambda        = 0.514;                          % wavelength in micron
ps            = 6.5/magnification;              % pixel size in micron
wavenumber    = 2*pi/lambda;                    % wave number
illu_rotation = [0, 180, 90, 270];              % orientation of the DPC half circular patterns           
num_rotation  = numel(illu_rotation);           % number of illumination used in DPC 
na_inner      = [0, 0, 0, 0];                   % if annular illumination is used, set the na corresponds to the inner radius     
SystemSetup();

%% show measurements
figure('Name', 'normalized, background substracted DPC measurements', 'NumberTitle', 'off')
for source_index = 1:num_rotation
   subplot(2, 2, source_index);
   imagesc(IDPC(:, :, source_index)); axis image; axis off; colormap gray; caxis([-0.5, 0.5]);
   title(['DPC ', num2str(source_index)], 'FontSize', 24);
end
drawnow;

%% generate illumination sources
%CHANGE THE SOURCE PATTERNS BASED ON YOUR ILLUMINATION SETTINGS
source      = zeros(dim(1), dim(2), num_rotation); 
for source_index = 1:num_rotation
    source(:, :, source_index) = genSource(illu_rotation(source_index), na_illum, na_inner(source_index), lambda, Fx, Fy);
end

figure('Name', 'LED Illumination Patterns', 'NumberTitle', 'off');
fig_rows    = floor(sqrt(num_rotation));
fig_cols    = floor((num_rotation)/fig_rows);
for fig_col = 1:num_rotation
    fig_index = fig_col;
    ax        = subplot(fig_rows, fig_cols, fig_index);
    imagesc(fftshift(fx), fftshift(fy), fftshift(source(:, :, fig_col)));
    axis image; axis off;
    title(['Source ', num2str(fig_col)]);
    colormap(ax, 'gray'); caxis([0, 0.1]);
end
drawnow;

%% DPC amplitude and phase recovery with Tikhonov regularization
% calculate frequency spectrum of the measurements
fIDPC           = F(IDPC);
% pupil function set by the numerical aperture of the imaging system
pupil           = (Fx.^2+Fy.^2<=(na/lambda)^2);
% parameters for Tikhonov regurlarization [absorption, phase] ((need to tune this based on SNR)
reg_tik         = 1.0*[1e-1, 5e-3];
% DPC deconvolution to solve the absorption and phase
[absorption,...
 phase]         = DPC_tik(fIDPC, source, pupil, reg_tik);

figure('Name', 'DPC Recovery Results', 'NumberTitle', 'off');
ax1 = subplot(1, 2, 1);
imagesc(x, y, absorption); axis image; axis off; colormap(ax1, 'gray'); caxis([-.15, 0.02]);
title('recovered \alpha','FontSize', 24);
ax2 = subplot(1, 2, 2);
imagesc(x, y, phase); axis image; axis off; colormap(ax2, 'gray'); caxis([-1.0, 3.0]);
title('recovered \phi', 'FontSize', 24);
linkaxes([ax1, ax2]);