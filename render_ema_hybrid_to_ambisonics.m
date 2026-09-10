% script uses discrete microphones, rigid-sphere HRTFs

clear;

addpath('dependencies/');

%input_file = 'wiener_filter/voice_recording_noisy.wav';
input_file = 'wiener_filter/voice_recording_denoised.wav';

output_file = 'out_ema_hybrid_binaural.wav';

% -------------------------------------------------------------------------

% array
N = 3; % microphone array
L = 7;   % no. of mic pairs
radial_filter_limit_dB = 50;

head_orientation_azimuth_deg = 0;

recorded_signals = audioread(input_file);

R_mic1 = 0.0565; % mic radius
R_mic2 = 0.0635; % mic radius

% filters
fc_int = 10;
order_int = 1;

% crossover between omni and virtual cardioid mic
fc_cross = 1000; % Hz

R_head = 0.08;

fs   = 48000;
taps = 2048*2;

c = 343;

s_mic_omni_1 = recorded_signals(:, 1:7); % inner ring of mics
s_mic_omni_2 = recorded_signals(:, 8:14); % outer ring of mics
% -------------------------------------------------------------------------

f     = linspace(0, fs/2, taps/2+1).';
omega = 2*pi*f;
k     = omega/c;
kr1   = k * R_mic1;
kr2   = k * R_mic2;

azi_mic = - (0:L-1) * 2*pi/L + pi; % radians

% ----------------- create cardioid signal from omnis ---------------------

% difference quotient
s_finite_diff = (s_mic_omni_2-s_mic_omni_1) / (R_mic2-R_mic1);
    
% choose omni signal
s_mic_omni = s_mic_omni_1;

[b, a] = get_highpassed_integrator(fs, fc_int, order_int, 'bilinear');

% signal of virtual cardioid
s_mic_card = s_mic_omni + c * filter(b, a, s_finite_diff);

% ------------------------ harmonic decomposition -------------------------
s_ring_nm_omni = get_sound_field_sh_coeffs_from_ema(s_mic_omni, N, azi_mic);
s_ring_nm_card = get_sound_field_sh_coeffs_from_ema(s_mic_card, N, azi_mic);

% --------------------------- get radial filters --------------------------
[b_m, b_nm_omni_t] = ema_2d_radial_filters(N, kr1, radial_filter_limit_dB, 'omni'); 
[ ~ , b_nm_card_t] = ema_2d_radial_filters(N, kr1, radial_filter_limit_dB, 'cardioid');

% ------------------------ Crossover radial filters -----------------------
[b_nm_omni_t, b_nm_card_t] = apply_crossover(b_nm_omni_t, b_nm_card_t, fs, fc_cross);

% -------------------------- Apply radial filters -------------------------
s_breve_nm_omni = zeros(size(b_nm_card_t, 1)+size(s_ring_nm_card, 1)-1, size(s_ring_nm_card, 2));
s_breve_nm_card = zeros(size(b_nm_card_t, 1)+size(s_ring_nm_card, 1)-1, size(s_ring_nm_card, 2));

for n = 0 : N 
    for m = -n : n        
        s_breve_nm_omni(:, n^2+n+m+1) = conv(b_nm_omni_t(:, n^2+n+m+1), s_ring_nm_omni(:, n^2+n+m+1));
        s_breve_nm_card(:, n^2+n+m+1) = conv(b_nm_card_t(:, n^2+n+m+1), s_ring_nm_card(:, n^2+n+m+1));
    end
end

% add crossover'd omni and cardioid coefficients to produce the final
% ambisonic signal
ambi_signal = s_breve_nm_omni + s_breve_nm_card;

% --------------------- Normalize the ambisonic signals -------------------

max_value   = max(abs(ambi_signal(:)));
weight      = 0.99 / max_value(1);
ambi_signal = ambi_signal .* weight;

% ---------------------- Store the ambisonic signals ----------------------

out_file_name = 'out_ambisonics.wav';  
%audiowrite(out_file_name, ambi_signal, fs);

% -------------------------- Binaural rendering ---------------------------

% SH coefficients of the HRIRs
hrirs = load(sprintf('resources/hrirs_ku100_magls_sh_N%d.mat', N));

binaural = render_ambi_signals_binaurally_custom_hrtf_t(ambi_signal, hrirs, head_orientation_azimuth_deg/180*pi, N);

binaural = 0.99 * binaural/max(abs(binaural(:)));

warning('Whoops, we gotta swap left and right channels of the binaural signal because I confused left and right during the recording.');
audiowrite(output_file, binaural(:, [2, 1]), fs);

fprintf('Binaural output signal written to file ''%s''\n.', output_file);




