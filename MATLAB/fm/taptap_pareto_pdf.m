function vPdf = taptap_pareto_pdf(vX,x_min,alpha1,x_c1,alpha2,x_c2,w)
%
%   Mixture of two tapered Pareto distributions
%
%   Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 18 June 2022
%   ...
%   version: 1.0.0, 18 June 2022
%
    vTap1 = tap_pareto_pdf(vX,alpha1,x_min,x_c1);
    vTap2 = tap_pareto_pdf(vX,alpha2,x_min,x_c2);
    vPdf = w*vTap1 + (1-w)*vTap2;
end
