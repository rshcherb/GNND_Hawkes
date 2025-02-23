function vPdf = trn_pareto_pdf(vX,gamma,x_min,x_max)
%
%   Truncated Pareto distribution: Kagan (GJI 2002 p.523)
%
%   Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 18 June 2022
%   ...
%   version: 1.0.0, 18 June 2022
%
    vPdf = vX.^(-1.0-gamma).*gamma./(x_min^(-gamma) - x_max^(-gamma));
end
