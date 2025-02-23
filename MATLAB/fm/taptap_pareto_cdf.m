function vPdf = taptap_pareto_cdf(vX,x_min,alpha1,x_c1,alpha2,x_c2,w)
%
%   CDF of the mixture of two tapered Pareto distributions
%
%   Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 18 November 2024
%   ...
%   version: 1.0.0, 18 November 2024
%
    vTap1 = tap_pareto_cdf(vX,alpha1,x_min,x_c1);
    vTap2 = tap_pareto_cdf(vX,alpha2,x_min,x_c2);
    vPdf = w*vTap1 + (1-w)*vTap2;
end
