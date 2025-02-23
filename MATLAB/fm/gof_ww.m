function [fGof, fRes] = gof_ww(vM,X_min,sFMD,vPar,OptArgs)
%
%   Computes the goodness-of-fit of a given distribution: Wiemer and Wyss, BSSA 90 (2000) p.859 
%
%   Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 23 June 2022
%   ...
%   version: 2.0.0, 18 November 2024
%
    arguments
        vM double                      % earthquake magnitudes or moments
        X_min double                   % lower magnitude or moment cutoff
        sFMD char                      % 'Exp', 'Pareto', 'TapPareto', 
        vPar double                    % parameters of the model distribution
        OptArgs.Units char = 'moment'  % in which units to plot: 'moment' or 'magnitude'
    end
    
    vM = vM(vM >= X_min);
    [cdf, x] = ecdf(vM);     % empirical cdf of the data
    B = 1 - cdf;             % complimentary cdf
    if strcmp(sFMD,'Exp')
        S = exp(-vPar(1)*(x - X_min)); % complimentary cdf and 
    elseif strcmp(sFMD,'Pareto')
        if strcmp(OptArgs.Units,'moment')
            S = 1 - pareto_cdf(x,vPar(1),X_min);
        elseif strcmp(OptArgs.Units,'magnitude')

        end
    elseif strcmp(sFMD,'TapPareto')
        if strcmp(OptArgs.Units,'moment')
            S = 1 - tap_pareto_cdf(x,vPar(1),X_min,vPar(2));
        elseif strcmp(OptArgs.Units,'magnitude')

        end
    elseif strcmp(sFMD,'TapTapPareto')
        if strcmp(OptArgs.Units,'moment')
            S = 1 - taptap_pareto_cdf(x,X_min,vPar(1),vPar(2),vPar(3),vPar(4),vPar(5));
        elseif strcmp(OptArgs.Units,'magnitude')

        end
    elseif strcmp(sFMD,'TapParetoPareto')
        if strcmp(OptArgs.Units,'moment')
            S = 1 - tap_pareto_pareto_cdf(x,X_min,vPar(1),vPar(2),vPar(3),vPar(4));
        elseif strcmp(OptArgs.Units,'magnitude')

        end
    end
    W = ones(length(B),1)/sum(B);    % weights
    fRes = sum(W.*abs(B - S));       % residual
    fGof = 100.0*(1.0 - fRes);       % goodness-of-fit
end
