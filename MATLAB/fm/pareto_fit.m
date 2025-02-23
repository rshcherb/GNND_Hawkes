function [vPar, vPci, vPerr] = pareto_fit(vData,M_min,OptArgs)
%
%   Fit the Pareto distribution to data: Kagan (GJI 2002 p.523)
%
%   Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 18 June 2022
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
    %                                                    beta
    if isempty(OptArgs.Start),      OptArgs.Start =      1.0; end
    if isempty(OptArgs.LowerBound), OptArgs.LowerBound = 0.0; end
    if isempty(OptArgs.UpperBound), OptArgs.UpperBound = 5.0; end

    options = statset('MaxIter',1000,'MaxFunEvals',500);
    [vPar, vPci] = mle(vData,'pdf',@(x,beta) pareto_pdf(x,beta,M_min),'Optimfun','fminsearch','Options',options,...
                      'LowerBound',OptArgs.LowerBound,'UpperBound',OptArgs.UpperBound,'Start',OptArgs.Start,'Alpha',OptArgs.Alpha);
    vPerr = 0.5*(vPci(2) - vPci(1));
end
