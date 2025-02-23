function [vPar_ci, vPar_err, mCovMat] = pp_param_ci(vPar_est,vCat,PP,OptArgs)
%
%   vCat  - the catalogue of events
%   PP    - input model structure: fT0, fTs, fTe, fT1, fM0, ...
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 31 October 2024
%   ...
%   version: 1.0.0, 2 November 2024
%
    arguments
        vPar_est double                       % the vector of the estimated PP parameters
        vCat double                           % event catalogue
        PP struct                             % structure that defines the poit-process to fit
        OptArgs.PointProcess char = 'Hawkes'  % 'ETAS'
        OptArgs.ErrorMethod char = 'Hessian'  % 'LogLik', 'BootStrap', 'Simul'
        %OptArgs.Hessian double = []           % provide the Hessian for teh model parameters
    end
    
    vPar_ci  = zeros(2,length(PP.mIP.vInitParams));
    vPar_err = zeros(1,length(PP.mIP.vInitParams));
    mCovMat  = zeros(length(vPar_est));
    % computes parameter errors for the PP model using one of the methods: 'Hesian', 'LogLik', etc
    if strcmp(OptArgs.ErrorMethod,'Simul')
        nSimul = 100;
        ci = ci_pp_simul(vPar_est,vCat,PP,nSimul,OptArgs.PointProcess);
        disp(ci)
    elseif strcmp(OptArgs.ErrorMethod,'Hessian') % using the numerical Hessian from the optimization
        mFI = PP.mHessian;
    elseif strcmp(OptArgs.ErrorMethod,'LogLik')
        mFI = [];
    elseif strcmp(OptArgs.ErrorMethod,'BootStrap')
    else
    end

    if strcmp(OptArgs.ErrorMethod,'BootStrap')
    else
        mCovMat  = inv(mFI); % inverse of the Fisher information matrix
        fPCalpha = norminv(1-0.5*PP.fAlpha,0,1); % (1-\alpha/2)*100% percentile of the standard normal distribution, N(0,1)
        for n = 1:PP.nPar
            vPar_err(n) = fPCalpha*sqrt(mCovMat(n,n));
        end
    end
end

