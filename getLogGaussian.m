function output = getLogGaussian( x, p )

output = -((p.q)*log(2*pi)) - (0.5*log(det(p.cov))) - 0.5*((x-p.m) * inv(p.cov) * (x-p.m)');
end

