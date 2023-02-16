clear;

% Converts a recording of an equatorial microphone array into ambisonic
% signals that can be rendered, for example with the tools from
%
%   https://plugins.iem.at/    
%   https://leomccormack.github.io/sparta-site/ .
%
% The original method was presented in 
%
%   J. Ahrens, H. Helmholz, D. L. Alon, S. V. Amengual Garí, "Spherical 
%   Harmonic Decomposition of a Sound Field Based on Observations Along the 
%   Equator of a Rigid Spherical Scatterer" in J. Acoust. Soc. Am. 150(2), 
%   2021, DOI: https://doi.org/10.1121/10.0005754.
%
% Detailed documentation of the present implementation is available in
% 
%   Jens Ahrens, "Ambisonic Encoding of Signals From Equatorial Microphone
%   Arrays," Technical note v.1, Chalmers University of Technology, 2022.
%   https://arxiv.org/abs/2211.00584
%
% (c) 2022 by Jens Ahrens

addpath('dependencies/');

% ----------------------------- Input data --------------------------------

% load the variables array_signals, fs, N, R, alpha_ema
load('resources/ema_recording_chalmers.mat');

% alpha_ema: azimuth of microphone positions in rad

% ----------------------------- Preparations ------------------------------

hankel_type = 2; 

radial_filter_length = 2048;

f = linspace(0, fs/2, radial_filter_length/2 + 1).'; 
c = 343;
k = 2*pi*f/c;

% -------------------- Precompute the radial filters ----------------------

gain_limit_radial_filters_dB = 40; % This is equivalent to 18 dB for the SMA.
reg_type_radial_filters = 'tikhonov';

[~, ema_inv_rf_t] = get_ema_radial_filters(k, R, N, gain_limit_radial_filters_dB, reg_type_radial_filters, hankel_type);

% ------------------------ Get the ambisonic signals ----------------------

ambi_signals = get_sound_field_sh_coeffs_from_ema_t(array_signals, ema_inv_rf_t, N, alpha_ema);

% --------------------- Normalize the ambisonic signals -------------------

max_value    = max(abs(ambi_signals(:)));
weight       = 0.99 / max_value(1);
ambi_signals = ambi_signals .* weight;

% ---------------------- Store the ambisonic signals ----------------------

out_file_name = 'out_ambisonics.wav';  
audiowrite(out_file_name, ambi_signals, fs);

% ------------------------ Create binaural preview ------------------------

% You would usually want to equalize the binaural rendering to mitigate 
% artifacts due to spherical harmonic order truncation and spatial
% aliasing. We'll add this in the future. 

head_orientation = 0;
out_lr = render_ambi_signals_binaurally_t(ambi_signals, head_orientation, N, 'transform_integral');

out_lr = out_lr / max(abs(out_lr(:)));
out_file_name = 'out_ema_binaural.wav';  
audiowrite(out_file_name, out_lr, fs);


