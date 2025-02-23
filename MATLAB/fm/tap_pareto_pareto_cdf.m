function vPdf = tap_pareto_pareto_cdf(vX,x_min,alpha,x_c,beta,w)
%
%   Mixture of two tapered Pareto distributions
%
%   Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 18 November 2024
%   ...
%   version: 1.0.0, 18 November 2024
%
    vTapPareto = tap_pareto_cdf(vX,alpha,x_min,x_c);
    vPareto = pareto_cdf(vX,beta,x_min);
    vPdf = w.*vTapPareto + (1-w).*vPareto;
end
