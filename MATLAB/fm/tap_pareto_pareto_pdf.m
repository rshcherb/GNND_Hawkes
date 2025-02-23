function vPdf = tap_pareto_pareto_pdf(vX,x_min,alpha,x_c,beta,w)
%
%   Mixture of two tapered Pareto distributions
%
%   Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 20 June 2022
%   ...
%   version: 1.0.0, 20 June 2022
%
    vTapPareto = tap_pareto_pdf(vX,alpha,x_min,x_c);
    vPareto = pareto_pdf(vX,beta,x_min);
    vPdf = w.*vTapPareto + (1-w).*vPareto;
end
