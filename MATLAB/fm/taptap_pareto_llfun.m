function fLLe = taptap_pareto_llfun(vMom,x_min,alpha1,x_c1,alpha2,x_c2,w)
%
%   Log-likelihood function for the mixture of the two tapered Pareto distributions
%
%   Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 16 November 2024
%   ...
%   version: 1.0.0, 16 November 2024
%
    vTap1 = tap_pareto_pdf(vMom,alpha1,x_min,x_c1);
    vTap2 = tap_pareto_pdf(vMom,alpha2,x_min,x_c2);
    vPdf = w*vTap1 + (1-w)*vTap2;

    fLLe = sum(log(vPdf));
end
