function [b_n, b_n_inv, b_n_inv_t] = get_sma_radial_filters(k, R, N, limit_db, reg_type, hankel_type)
% b_n is given by Eq. (6). b_n_inv and b_n_inv_t are the inverse of that in
% frequency domain and time domain, respectively.
%
% reg_type: 'tikhonov', 'hard', 'soft'
%
% b_n_inv and b_n_inv_t are causal (they comprise a modeling delay).
%
% The equation numbers refer to 
% 
%   Jens Ahrens, "Ambisonic Encoding of Signals From Spherical Microphone
%   Arrays," Technical note v. 1, Chalmers University of Technology, 2022.
%   https://arxiv.org/abs/2211.00583
%
% Written by Jens Ahrens, 2022

% ----------------------- compute the terms b_n, Eq. (6) ------------------
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

   for index = indices
       % this usually occurs at the DC bin, so replace that one.
       b_n(index) = abs(b_n(index+1));    
   end
   
end

% invert
b_n_inv = 1./b_n;

% ----------------------------- limiting ----------------------------------
if limit_db < inf
    
    % --- soft clipping ---
    if strcmp(reg_type, 'soft')
    
        b_n_inv = (2*10^(limit_db/20)) / pi * b_n_inv ./ abs(b_n_inv) .* atan(pi/(2*10^(limit_db/20)) .* abs(b_n_inv));
        
    % --- hard clipping ---
    elseif strcmp(reg_type, 'hard')
        
        b_n_inv(abs(b_n_inv) > 10^(limit_db/20)) = b_n_inv(abs(b_n_inv) > 10^(limit_db/20)) ./ abs(b_n_inv(abs(b_n_inv) > 10^(limit_db/20))) * 10^(limit_db/20); 
        
    % --- Tikhonov regularization as used in (Moreau, Daniel, Bertet, AES 2006) ---
    elseif strcmp(reg_type, 'tikhonov')
        
       limit          = 10^(limit_db/20);
       lambda_squared = (1 - sqrt(1 - 1/limit^2)) ./ (1 + sqrt(1 - 1/limit^2));
       
       b_n_inv     = conj(b_n) ./ (abs(b_n).^2 + lambda_squared);      
    
    else
        
        error('Unknown reg_type.');
        
    end

end

% ----- catch possible extra NaNs that can occur during the inversion -----
indices = find(isnan(b_n_inv));

if ~isempty(indices)
    warning('%d more occurrence(s) of NaN deteced in high-order radial filters. Applying a fix.', length(indices));

   for index = indices
       % this usually occurs at the DC bin, so replace that one.
       b_n_inv(index) = abs(b_n_inv(index+1));    
   end
   
end

% --------------------------- make filters causal -------------------------

b_n_inv_sym = [b_n_inv; conj(flipud(b_n_inv(2:end-1, :)))];

b_n_inv_t = ifft(b_n_inv_sym, 'symmetric'); 

filter_length = size(b_n_inv_t, 1);

% make causal (we assume that the filter_length is 2048 or longer)
b_n_inv_t = circshift(b_n_inv_t, [filter_length/2 0]);

b_n_inv = fft(b_n_inv_t);
b_n_inv = b_n_inv(1:end/2+1, :);

end


