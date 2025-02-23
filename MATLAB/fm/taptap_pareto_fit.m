function [vPar, vPci, vPerr] = taptap_pareto_fit(vData,M_min,OptArgs)
%
%   Fit the Taper. Pareto + Taper. Pareto distribution to the data: Kagan (GJI 2002 p.523)
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
    %                                                    [alpha1, Mc_1,       alpha2, Mc_2,      w]
    if isempty(OptArgs.Start),      OptArgs.Start =      [0.7     1.01*M_min  0.5     4*M_min    0.5]; end
    if isempty(OptArgs.LowerBound), OptArgs.LowerBound = [0.0     M_min       0.0     1.1*M_min  0.0]; end
    if isempty(OptArgs.UpperBound), OptArgs.UpperBound = [5.0     1e26        5.0     1e26       1.0]; end

    pdffun = @(x,alpha1,x_c1,alpha2,x_c2,w) taptap_pareto_pdf(x,M_min,alpha1,x_c1,alpha2,x_c2,w);
    %logpdffun = @(x,alpha1,x_c1,alpha2,x_c2,w) taptap_pareto_logpdf(x,x_min,alpha1,x_c1,alpha2,x_c2,w);
    
    options = statset('MaxIter',1000,'MaxFunEvals',2000);
    %options = statset('MaxIter',1000,'MaxFunEvals',2000,'Display','iter');
    [vPar, vPci] = mle(vData,'pdf',pdffun,'OptimFun','fminsearch','Options',options,...
                     'LowerBound',OptArgs.LowerBound,'UpperBound',OptArgs.UpperBound,'Start',OptArgs.Start,'Alpha',OptArgs.Alpha);
%     [vPar, vPci] = mle(vData,'logpdf',logpdffun,'Options',options,...
%                      'LowerBound',OptArgs.LowerBound,'UpperBound',OptArgs.UpperBound,'Start',OptArgs.Start,'Alpha',OptArgs.Alpha);

    vPerr(1) = 0.5*(vPci(2,1) - vPci(1,1));
    vPerr(2) = 0.5*(vPci(2,2) - vPci(1,2));
    vPerr(3) = 0.5*(vPci(2,3) - vPci(1,3));
    vPerr(4) = 0.5*(vPci(2,4) - vPci(1,4));
    vPerr(5) = 0.5*(vPci(2,5) - vPci(1,5));
end
