function vPdf = taptap_pareto_logpdf(vX,x_min,alpha1,x_c1,alpha2,x_c2,w)
%
%   Mixture of two tapered Pareto distributions
%
%   Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 18 June 2022
%   ...
%   version: 1.0.0, 18 June 2022
%
%     logTap1 = tap_pareto_logpdf(vX,alpha1,x_min,x_c1);
%     logTap2 = tap_pareto_logpdf(vX,alpha2,x_min,x_c2);
%     %vPdf = w*vTap1 + (1-w)*vTap2;
%     vPdf = log(w) + logTap1 + log(1-w) + logTap2;
    vPdf = log(w*(x_min./vX).^alpha1.*(alpha1./vX + 1.0/x_c1).*exp((x_min - vX)/x_c1) + ...
               (1.0-w)*(x_min./vX).^alpha2.*(alpha2./vX + 1.0/x_c2).*exp((x_min - vX)/x_c2));
end
