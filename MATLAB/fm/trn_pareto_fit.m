function [vPar, vPci, vPerr] = trn_pareto_fit(vData,M_min,OptArgs)
%
%   Fit the truncated Pareto distribution to data: Kagan (GJI 2002 p.523)
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
    %                                                    [gamma, M_max]
    if isempty(OptArgs.Start),      OptArgs.Start =      [1.1    10*M_min]; end
    if isempty(OptArgs.LowerBound), OptArgs.LowerBound = [0.01   5*M_min]; end
    if isempty(OptArgs.UpperBound), OptArgs.UpperBound = [10.0   1e26]; end

    options = statset('MaxIter',1000,'MaxFunEvals',500);
    [vPar, vPci] = mle(vData,'pdf',@(x,gamma,x_max) trn_pareto_pdf(x,gamma,M_min,x_max),'Optimfun','fminsearch','Options',options,...
                    'LowerBound',OptArgs.LowerBound,'UpperBound',OptArgs.UpperBound,'Start',OptArgs.Start,'Alpha',OptArgs.Alpha);
    vPerr(1) = 0.5*(vPci(2,1) - vPci(1,1));
    vPerr(2) = 0.5*(vPci(2,2) - vPci(1,2));
end
