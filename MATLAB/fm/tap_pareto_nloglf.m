function nloglf = tap_pareto_nloglf(vP,vData,x_min)
%
%   Negative log-likelihood of the tapered Pareto distribution: Kagan (GJI 2002 p.525)
%
%   Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 18 June 2022
%   ...
%   version: 1.0.0, 18 June 2022
%
    alpha = vP(1);
    x_c = vP(2);
    n = length(vData);
    nloglf = -(n.*alpha.*log(x_min) + (n.*x_min - sum(vData))./x_c - alpha.*sum(log(vData)) + sum(log(alpha./vData + 1.0./x_c)));
    %disp(nloglf)
end
