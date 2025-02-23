function vPdf = pareto_pdf(vX,beta,x_min)
%
%   Tapered Pareto distribution: Kagan (GJI 2002 p.521)
%
%   Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 18 June 2022
%   ...
%   version: 1.0.0, 18 June 2022
%
    vPdf = beta*x_min^beta*vX.^(-1.0-beta); % Kagan (GJI 2002 p.521), Eq. (4), Zaliapin et al (PAGEOPH 2005 p.1189)
end
