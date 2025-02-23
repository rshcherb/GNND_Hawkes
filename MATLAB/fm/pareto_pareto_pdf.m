function vPdf = pareto_pareto_pdf(vX,x_min,gamma,x_c,beta)
%
%   Mixture of two tapered Pareto distributions
%
%   Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 21 June 2022
%   ...
%   version: 1.0.0, 21 June 2022
%
    indx = vX >= x_c;
    vPdf = zeros(size(vX));
    A = beta*gamma/(gamma*x_c^(-beta) + beta*(x_min^(-gamma) - x_c^(-gamma)));

    vPdf(indx) = A.*vX(indx).^(-1.0-beta); % Pareto above x_c
    vPdf(~indx) = A.*vX(~indx).^(-1.0-gamma); % truncated Pareto between x_min and x_c
end
