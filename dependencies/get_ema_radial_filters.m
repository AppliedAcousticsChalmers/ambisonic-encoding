function [ema_inv_rf, ema_inv_rf_t] = get_ema_radial_filters(k, R, N, limit_db, reg_type, hankel_type)
% ema_inv_rf and ema_inv_rf_t are the 1st factor on the right hand side of
% Eq. (13) in frequency and time domain respectively. 
%
% reg_type: 'tikhonov', 'hard', 'soft'
%
% They are causal (they comprise a modeling delay).
%
% The equation numbers refer to 
% 
%   Jens Ahrens, "Ambisonic Encoding of Signals From Equatorial Microphone
%   Arrays," Technical note v. 1, Chalmers University of Technology, 2022.
%   https://arxiv.org/abs/2211.00584
%
% Written by Jens Ahrens, 2022

% ----------------------- Compute the terms b_n, Eq. (7) ------------------
b_n = zeros(size(k, 1), N+1, 'like', 1j);

kR = k.*R + 5*eps; % add 5*eps to avoid undefined values for the Hankel function

for n = 0 : N

    hankel_prime = 1/(2*n+1) * (n * sphbesselh(n-1, hankel_type, kR) - (n+1) * sphbesselh(n+1, hankel_type, kR));

    b_n(:, n+1) = 1i./kR.^2 .* 1./hankel_prime;

    if hankel_type == 1
        b_n(:, n+1) =   4*pi * 1i^(-n) .* b_n(:, n+1);
    elseif hankel_type == 2
        b_n(:, n+1) = - 4*pi * 1i^n    .* b_n(:, n+1);
    else 
        error('Unknown hankel_type.');
    end

end

% ----------------- catch possible NaNs (can occur for N > 19) ------------
indices = find(isnan(b_n));

if ~isempty(indices)

    warning('%d occurrence(s) of NaN deteced in high-order radial filters. Applying a fix.', length(indices));
    % this usually occurs at the DC bin, so replace that one.
    b_n(indices) = abs(b_n(indices+1));

end

% ------------ Compute the inverse EMA radial filters, Eq. (13) -----------
ema_rf = zeros(size(b_n, 1), 2*N+1, 'like', b_n);

for m = -N : N

    for n_prime = abs(m) : N

        ema_rf(:, m+N+1) = ema_rf(:, m+N+1) + b_n(:, n_prime+1) .* N_nm(n_prime, m, pi/2).^2;

    end

end

% -------------------------------invert -----------------------------------
ema_inv_rf = 1./ema_rf;

% ----------------------------- limiting ----------------------------------
if limit_db < inf

    limit = 10^(limit_db/20);

    % --- soft clipping ---
    if strcmpi(reg_type, 'soft')

         ema_inv_rf = (2*limit)/pi * ema_inv_rf ./ abs(ema_inv_rf) .* atan(pi/(2*limit) .* abs(ema_inv_rf));

    % --- hard clipping ---
    elseif strcmpi(reg_type, 'hard')

        indices = abs(ema_inv_rf) > limit;
        ema_inv_rf(indices) = ema_inv_rf(indices) ./ abs(ema_inv_rf(indices)) * limit;

    % --- Tikhonov regularization as used in (Moreau, Daniel, Bertet, AES 2006) ---
    elseif strcmpi(reg_type, 'tikhonov')

        lambda_squared = (1 - sqrt(1 - 1/limit^2)) ./ (1 + sqrt(1 - 1/limit^2));
        ema_inv_rf = conj(ema_rf) ./ (abs(ema_rf).^2 + lambda_squared);

    else

        error('Unknown reg_type.');

    end

end

% ----- catch possible extra NaNs that can occur during the inversion -----
indices = find(isnan(ema_inv_rf));

if ~isempty(indices)

    warning('%d more occurrence(s) of NaN deteced in high-order radial filters. Applying a fix.', length(indices));
    % this usually occurs at the DC bin, so replace that one.
    ema_inv_rf(indices) = abs(ema_inv_rf(indices+1));

end

% --------------------------- make filters causal -------------------------
ema_inv_rf_sym = [ema_inv_rf; conj(flipud(ema_inv_rf(2:end-1, :)))];
ema_inv_rf_t = ifft(ema_inv_rf_sym, 'symmetric');

% make causal (we assume that the filter_length is 2048 or longer)
ema_inv_rf_t = circshift(ema_inv_rf_t, size(ema_inv_rf_t, 1)/2);

ema_inv_rf = fft(ema_inv_rf_t);
ema_inv_rf = ema_inv_rf(1:end/2+1, :);

end
