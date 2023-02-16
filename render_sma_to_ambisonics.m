clear;

% Converts a recording of a rigid spherical microphone array into ambisonic
% signals that can be rendered, for example with the tools from
%
%   https://plugins.iem.at/    
%   https://leomccormack.github.io/sparta-site/ .
%
% A detailed documentation is available in 
% 
%   Jens Ahrens, "Ambisonic Encoding of Signals From Spherical Microphone
%   Arrays," Technical note v. 1, Chalmers University of Technology, 2022.
%   https://arxiv.org/abs/2211.00583
%
% (c) 2022 by Jens Ahrens

addpath('dependencies/');

% ----------------------------- Input data --------------------------------

% load the variables array_signals, fs, N, R, alpha_sma, beta_sma
load('resources/eigenmike_walking_around.mat');

% alpha_sma: azimuth of microphone positions in rad
% beta_sma: colatitude (= polar angle) of microphone positions in rad

% ----------------------------- Preparations ------------------------------

hankel_type = 2; 

radial_filter_length = 2048;

f = linspace(0, fs/2, radial_filter_length/2 + 1).'; 
c = 343;
k = 2*pi*f/c;

% ---------------------- Precompute radial filters ------------------------

gain_limit_radial_filters_dB = 20;
reg_type_radial_filters = 'tikhonov';

[~, ~, b_n_inv_t] = get_sma_radial_filters(k, R, N, gain_limit_radial_filters_dB, reg_type_radial_filters, hankel_type);

% ------------------------ Get the ambisonic signals ----------------------

ambi_signals = get_sound_field_sh_coeffs_from_sma_t(array_signals, b_n_inv_t, N, beta_sma, alpha_sma, grid_weights_array);

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

out_lr = out_lr / max(abs(out_lr(:))) * 0.99; % normalize the binaural signal
out_file_name = 'out_sma_binaural.wav';  
audiowrite(out_file_name, out_lr, fs);



