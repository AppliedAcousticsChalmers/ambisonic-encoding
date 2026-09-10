function [out] = besselh_prime(m, nu, x)

out = 0.5 * (besselh(m-1, nu, x) - besselh(m+1, nu, x));

end

