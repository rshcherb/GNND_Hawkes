function [vPar, vPci, vPerr] = tap_pareto_fit(vData,M_min,OptArgs)
%
%   Fit the tapered Pareto distribution to data: Kagan (GJI 2002 p.523)
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
    %                                                    [alpha,  x_c0]
    if isempty(OptArgs.Start),      OptArgs.Start =      [1.4     11*M_min]; end
    if isempty(OptArgs.LowerBound), OptArgs.LowerBound = [0.0     M_min]; end
    if isempty(OptArgs.UpperBound), OptArgs.UpperBound = [10.0    1e26]; end

    %pdffun = @(x,alpha,xc) tap_pareto_pdf(x,alpha,M_min,xc);
    %cdffun = @(x,alpha,xc) tap_pareto_cdf(x,alpha,M_min,xc);
    %logpdffun = @(x,alpha,xc) tap_pareto_logpdf(x,alpha,M_min,xc);
    nloglfun = @(vP,vData,cens,freq) tap_pareto_nloglf(vP,vData,M_min);

    options = statset('MaxIter',1000,'MaxFunEvals',500);
%     try
%         [vPar, vPci] = mle(vData,'pdf',pdffun,'Options',options,'LowerBound',OptArgs.LowerBound,'Start',OptArgs.Start,'Alpha',OptArgs.Alpha);
%         [vPar, vPci] = mle(vData,'pdf',pdffun,'cdf',cdffun,'TruncationBounds',[M_min Inf],'Options',options,...
%                          'LowerBound',OptArgs.LowerBound,'Start',OptArgs.Start,'Alpha',OptArgs.Alpha);
%         [vPar, vPci] = mle(vData,'logpdf',logpdffun,'Options',options,'LowerBound',OptArgs.LowerBound,'Start',OptArgs.Start,'Alpha',OptArgs.Alpha);
        [vPar, vPci] = mle(vData,'nloglf',nloglfun,'Optimfun','fminsearch','Options',options,...
            'LowerBound',OptArgs.LowerBound,'UpperBound',OptArgs.UpperBound,'Start',OptArgs.Start,'Alpha',OptArgs.Alpha);
%     catch ME
%         disp(ME.message)
%     end

    % using Kagan's estimates
%     [vPar, vPci] = vere_jones_mle(vData,X_min);

    vPerr(1) = 0.5*(vPci(2,1) - vPci(1,1));
    vPerr(2) = 0.5*(vPci(2,2) - vPci(1,2));
end

function [vPar, vPci] = vere_jones_mle(vData,M_min)
%
% Vere-Jones et al. (2001), Kagan (GJI 2002 p.541)
%
    n = length(vData);
    S = vData/M_min;
    A = 1.0/n*sum(log(S)); % Eq. (A6)
    B = 1.0/n*sum(S-1.0);  % Eq. (A7)
    eta2 = 1.0/(1.01*M_min);
    %eta2 = 1.0/(B - A*min(S));
    disp(eta2)
    eta1 = 0;
    myfun = @(x,S,n,A,B) (sum(1.0./(1.0 - x.*(B - A.*S))) - double(n));
    fun = @(x) myfun(x,S,n,A,B);
    eta_hat = fzero(fun,[eta1, eta2]);
    alpha_hat = (1.0 - eta_hat*B)/A;

    vPar(1) = alpha_hat;
    vPar(2) = 1/eta_hat;
    vPci = zeros(2,2);
end

