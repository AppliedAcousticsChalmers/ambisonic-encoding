clear;

% Low value: less aggressive noise suppression; large value: more aggressive noise suppression
diagonal_loading = 1e-4; 

% --- raw microphone array signals ---
noise_file  = 'background_noise_10s.wav';
input_file  = 'voice_recording_noisy.wav';

% --- denoised microphone array signals --- 
output_file = 'voice_recording_denoised.wav';

win_length = 1024;
hop_size   = 512;
nfft       = win_length;

% ---------------------------- load signals -------------------------------

% noisy signal
[sig_in, fs] = audioread(input_file);
[~, no_of_ch] = size(sig_in);

fprintf('Using noise signature from file ''%s''.\n', noise_file);

% noise
[noise, fs_noise] = audioread(noise_file);

if fs_noise ~= fs
    error('Noise file sampling rate does not match noisy file.');
end

if size(noise,2) ~= no_of_ch
    error('Noise file must have the same number of channels.');
end


% -------------------------------- STFT -----------------------------------

win = sqrt(hann(win_length,'periodic'));
no_of_bins = nfft/2 + 1;

[sig_spec,   no_of_sig_frames]   = stft(sig_in, win, hop_size, nfft);
[noise_spec, no_of_noise_frames] = stft(noise , win, hop_size, nfft);

% -------------------- estimate covariance matrices -----------------------

Rxx = zeros(no_of_ch, no_of_ch, no_of_bins);
Rnn = zeros(no_of_ch, no_of_ch, no_of_bins);

for bin = 1 : no_of_bins

    Xf = squeeze(sig_spec(  bin, :, :)).';   % channels x frames
    Nf = squeeze(noise_spec(bin, :, :)).';   % channels x frames

    Rxx(:, :, bin) = (Xf * Xf') / no_of_sig_frames;
    Rnn(:, :, bin) = (Nf * Nf') / no_of_noise_frames;

end

% ---------------------- multichannel Wiener filter -----------------------

out_spec = zeros(no_of_bins,no_of_sig_frames,no_of_ch);

for bin = 1 : no_of_bins

    Rx = Rxx(:, :, bin);
    Rn = Rnn(:, :, bin);

    Rx = 0.5 * (Rx + Rx');
    Rn = 0.5 * (Rn + Rn');

    Rs = Rx - Rn;
    Rs = 0.5 * (Rs + Rs');

    % Force Rs to be positive semidefinite
    [V, D] = eig(Rs);
    d = real(diag(D));
    d(d < 0) = 0;
    Rs = V * diag(d) * V';

    % Diagonal loading
    loading = diagonal_loading * trace(Rx)/no_of_ch;
    Rx_loaded = Rx + loading * eye(no_of_ch);

    % Wiener filter matrix
    W = Rs / Rx_loaded;

    % apply filter
    for frame = 1 : no_of_sig_frames
        xft = squeeze(sig_spec(bin, frame, :));  % channels x 1
        yft = W * xft;                         % channels x 1

        out_spec(bin, frame, :) = yft;
    end
end

% ----------------------------- inverse STFT ------------------------------

sig_out = istft(out_spec, win, hop_size, nfft);

% ----------------------------- store output ------------------------------

% TODO: fix first frame
%sig_out = sig_out(win_length+1:end, :);
 
audiowrite(output_file, sig_out, fs, 'BitsPerSample', 32);

fprintf('Denoised array signals written to ''%s''.\n', output_file);



% -------------------------------------------------------------------------
% --------------------------- Local functions -----------------------------
% -------------------------------------------------------------------------

function [X, no_of_frames] = stft(x, win, hop_size, nfft)

    [no_of_samples, no_of_channels] = size(x);
    win_length = length(win);

    no_of_frames = floor((no_of_samples - win_length)/hop_size) + 1;
    no_of_bins = nfft/2 + 1;

    X = zeros(no_of_bins, no_of_frames, no_of_channels);

    for ch = 1 : no_of_channels
        for frame = 1 : no_of_frames
            
            idx = (1:win_length) + (frame-1)*hop_size;

            x_frame = x(idx,ch) .* win;
            Xtmp = fft(x_frame, nfft);

            X(:, frame, ch) = Xtmp(1:no_of_bins);

        end
    end
end

function y = istft(Y, win, hop_size, nfft)

    [~, no_of_frames, no_of_channels] = size(Y);
    win_length = length(win);

    no_of_samples = (no_of_frames-1) * hop_size + win_length;
    y = zeros(no_of_samples, no_of_channels);

    for ch = 1 : no_of_channels

        % current channel
        y_ch = zeros(no_of_samples, 1);

        for frame = 1 : no_of_frames

            idx = (1:win_length) + (frame-1) * hop_size;

            Ypos = Y(:,frame,ch);

            Yfull = [Ypos; conj(Ypos(end-1:-1:2))];

            y_frame = real(ifft(Yfull,nfft));
            y_frame = y_frame(1:win_length) .* win;

            % add frame to output buffer
            y_ch(idx) = y_ch(idx) + y_frame;

        end
    
        % add channel to output
        y(:,ch) = y_ch;

    end

end