function [out] = besselj_prime(m, x)

out = 0.5 * (besselj(m-1, x) - besselj(m+1, x));

end

