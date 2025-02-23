function vCdf = pareto_cdf(vX,beta,x_min)
%
%   CDF of the Pareto distribution: Kagan (GJI 2002 p.521)
%
%   Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 18 November 2024
%   ...
%   version: 1.0.0, 18 November 2024
%
    vCdf = 1.0 - (x_min./vX).^beta; % Kagan (GJI 2002 p.523), Eq. (5), Zaliapin et al (PAGEOPH 2005 p.1189)
end
