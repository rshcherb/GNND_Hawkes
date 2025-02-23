function PP = model_etas_rate(vCat,PP,OptArgs)
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 23 May, 2017
%   ...
%   version: 4.0.0, 3 November 2024
%
    arguments
        vCat double
        PP struct
        OptArgs.PointProcess char    = 'ETAS'       % 'ETASFrack', 'ETASFrackConv', 'ETAS_RS', 
        OptArgs.DateFormat char      = 'days'       % to use for time axis as 'days' or 'date'
        OptArgs.ErrorMethod char     = 'Hessian'    % 'LogLik', 'BootStrap', 'Simul'
        OptArgs.OptimMethod char     = 'fmincon'    % 'gs'; 'fminsearch';
        OptArgs.BinNumber double     = []           % number of binns
        OptArgs.PlotRate char        = 'on'         % 
        OptArgs.PlotCumulative char  = 'on'         % 
        OptArgs.PlotTransformed char = 'on'         % 
        OptArgs.PlotCumulSum char    = 'on'         % 
        OptArgs.Title char           = []           % 
        OptArgs.SaveFigure char      = []           % 
        OptArgs.Position double = [400,300,800,600] % figure size and position
    end

    indx = (vCat(:,1) >= PP.fT0) & (vCat(:,1) <= PP.fTe);
    vCat = vCat(indx,:);
    if strcmp(OptArgs.PointProcess,'ETASFrack') || strcmp(OptArgs.PointProcess,'ETASFrackConv') || strcmp(OptArgs.PointProcess,'ETAS_MultConv') ||...
            strcmp(OptArgs.PointProcess,'ETAS_RS')
        indx = (PP.vT >= PP.fT0) & (PP.vT <= PP.fTe);
        PP.vT = PP.vT(indx);
    end
    vCat(:,2) = vCat(:,2) - PP.fM0; % subtract the reference magnitude fM0
    
    if strcmp(OptArgs.PointProcess,'ETAS')
        gif_pp = @(t,vCat,vPar,bRateType) gif_etas(t,vCat,PP.fTs,PP.nJs,vPar,bRateType);
        negLLfun = @(vPar) -llfun_etas(vPar,vCat,PP.fTs,PP.fTe,PP.nJs,PP.nJe);
    elseif strcmp(OptArgs.PointProcess,'ETASFrack')
        negLLfun = @(vPar) -llfun_etas_frack(vPar,vCat,PP.fTs,PP.fTe,PP.nJs,PP.nJe,PP.vT,PP.vTInjRate);
    elseif strcmp(OptArgs.PointProcess,'ETASFrackConv')
        gif_pp = @(t,vCat,vPar,bRateType) gif_etas_frack_conv(t,vCat,PP.vTInjRate,PP.fTs,PP.nJs,vPar,PP.nConvKernel,bRateType);
        negLLfun = @(vPar) -llfun_etas_frack_conv(vPar,vCat,PP.fTs,PP.fTe,PP.nJs,PP.nJe,PP.vTInjRate,PP.nConvKernel);
    elseif strcmp(OptArgs.PointProcess,'ETAS_MultConv')
        negLLfun = @(vPar) -llfun_etas_multconv(vPar,vCat,PP.fTs,PP.fTe,PP.nJs,PP.nJe,PP.vTInjRate,PP.nConvKernel);
    elseif strcmp(OptArgs.PointProcess,'ETAS_RS')
        negLLfun = @(vPar) -llfun_etas_rs(vPar,vCat,PP.fT0,PP.fTs,PP.fTe,PP.nJs,PP.nJe,PP.vT,PP.vTInjRate);
    end

    tic
    [PP.vPar, PP.fLle, exitflag, PP.mHessian] = pp_fit_mle(negLLfun,PP,'OptimMethod',OptArgs.OptimMethod);
    PP.nPar = length(PP.vPar);
    PP.AIC = -2.0*PP.fLle + 2*PP.nPar; % standard definition
    [PP.vPar_ci, PP.vPar_err, PP.mCovMat] = pp_param_ci(PP.vPar,vCat,PP);
    toc
    disp(['Exit flag: ',num2str(exitflag)]);

    if ~isempty(OptArgs.SaveFigure)
        save([OptArgs.SaveFigure,'_eq_cat.dat'],'vCat','-ascii');
        vParSave = [PP.vPar; PP.vPar_err];
        save([OptArgs.SaveFigure,'_parameters.txt'],'vParSave','-ascii');
    end
    
    plot_pp_rate(vCat,PP,gif_pp,...
        'PlotRate',OptArgs.PlotRate,'PlotCumulative',OptArgs.PlotCumulative,'PlotTransformed',OptArgs.PlotTransformed,'PlotCumulSum',OptArgs.PlotCumulSum,...
        'BinNumber',OptArgs.BinNumber,'DateFormat',OptArgs.DateFormat,'Title',OptArgs.Title,'SaveFigure',OptArgs.SaveFigure,'Position',OptArgs.Position);
end

