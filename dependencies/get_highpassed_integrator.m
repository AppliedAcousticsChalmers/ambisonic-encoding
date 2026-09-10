function [b, a] = get_highpassed_integrator(fs, fc, order, method)
% Integrator combined with a Butterworth highpass 
% order: order of highpass
% Method: 'bilinear' or 'impinv'

omega_c = 2*pi*fc;

% --------------------------- arbitrary order -----------------------------

% Analog Butterworth highpass
[num_hp, den_hp] = butter(order, omega_c, 'high', 's');

% Combine with ideal integrator 1/s

% H(s) = H_hp(s) / s
num_s = num_hp;
den_s = conv(den_hp, [1 0]);   % multiply denominator by s

% ---------------- Convert analog filter to digital filter ----------------

% bilinear is preferrable
if strcmp(method, 'bilinear')
    [b, a] = bilinear(num_s, den_s, fs);
elseif strcmp(method, 'impinv') 
    [b, a] = impinvar(num_s, den_s, fs);
else
    error('Unknown method.')
end

end