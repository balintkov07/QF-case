function logpdf = Logpdf_normal(y, mu, h)
% Logpdf_normal evaluates the log-density of a normal distribution.

logpdf = -0.5 * log(2*pi) - 0.5 * log(h) - 0.5 * (y - mu)^2 / h;

end
