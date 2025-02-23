function [vPar, vPci, vPerr] = pareto_pareto_fit(vData,M_min,OptArgs)
%
%   Fit the Pareto distribution to data: Kagan (GJI 2002 p.523)
%
%   Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 21 June 2022
%   ...
%   version: 2.0.0, 17 November 2024
%
    arguments
        vData double                   % data to fit the distribution
        M_min double                   % lower bound for the distribution values
        OptArgs.Start double = []      % initial parameters for the model
        OptArgs.LowerBound double = [] % lower bounds for the model parameters
        OptArgs.UpperBound double = [] % upper bounds for the model parameters
        OptArgs.Alpha double           % confidence level
    end
    %                                                    [gamma,  Mc,         beta]
    if isempty(OptArgs.Start),      OptArgs.Start =      [0.5     100*M_min   1.4]; end
    if isempty(OptArgs.LowerBound), OptArgs.LowerBound = [0.0     2*M_min     1.0]; end
    if isempty(OptArgs.UpperBound), OptArgs.UpperBound = [5.0     1e26        5.0]; end

    pdffun = @(x,gamma,x_c,beta) pareto_pareto_pdf(x,M_min,gamma,x_c,beta);
    
    options = statset('MaxIter',1000,'MaxFunEvals',1000);
    %options = statset('MaxIter',1000,'MaxFunEvals',2000,'Display','iter');
    [vPar, vPci] = mle(vData,'pdf',pdffun,'Optimfun','fminsearch','Options',options,...
                     'LowerBound',OptArgs.LowerBound,'UpperBound',OptArgs.UpperBound,'Start',OptArgs.Start,'Alpha',OptArgs.Alpha);

    vPerr(1) = 0.5*(vPci(2,1) - vPci(1,1));
    vPerr(2) = 0.5*(vPci(2,2) - vPci(1,2));
    vPerr(3) = 0.5*(vPci(2,3) - vPci(1,3));
end
