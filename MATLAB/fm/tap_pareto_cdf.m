function vCdf = tap_pareto_cdf(vX,alpha,x_min,x_c)
%
%   Tapered Pareto distribution: Kagan (GJI 2002 p.523)
%
%   Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 19 June 2022
%   ...
%   version: 1.0.0, 19 June 2022
%
    vCdf = 1.0 - (x_min./vX).^alpha.*exp((x_min - vX)/x_c);
end
