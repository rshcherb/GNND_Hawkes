function vPdf = tap_pareto_logpdf(vX,alpha,x_min,x_c)
%
%   Tapered Pareto distribution: Kagan (GJI 2002 p.523)
%
%   Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 18 June 2022
%   ...
%   version: 1.0.0, 18 June 2022
%
    %vPdf = (x_min./vX).^alpha.*(alpha./vX + 1.0/x_c).*exp((x_min - vX)/x_c);
    vPdf = alpha*log(x_min./vX) + log(alpha./vX + 1.0/x_c) + (x_min - vX)/x_c;
end
