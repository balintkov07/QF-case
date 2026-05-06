function [NegativeLogLikelihood, h] = GARCH_Filter(parameter_vector, returns)
% GARCH_Filter applies the GARCH(1,1) recursion and computes the
% negative Gaussian log-likelihood.
%
% parameter_vector = [mu; omega; alpha; beta]

T = length(returns);

mu    = parameter_vector(1);
omega = parameter_vector(2);
alpha = parameter_vector(3);
beta  = parameter_vector(4);


% The unconditional variance is only well-defined when alpha + beta < 1.
% If this condition is violated, return an infinite objective value.
if alpha + beta >= 1
    NegativeLogLikelihood = 10^10; % set to a very large value if the constraint is violated
    h = NaN(T,1);
    return
end

% Prefill
loglikelihood = zeros(T,1);
h = zeros(T,1);

% Initialise the conditional variance.
h(1) = omega / (1 - alpha - beta);

% Loop over time
for t = 1:T
   loglikelihood(t) = Logpdf_normal(returns(t), mu, h(t));
   if t<T
    h(t+1) = omega + alpha * (returns(t) - mu)^2 + beta * h(t);
   end
end

% Compute log likelihood
NegativeLogLikelihood = -sum(loglikelihood(1:end));

end
