function PP = model_hawkes_rate(vCat,PP,OptArgs)
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 31 October 2024
%   ...
%   version: 1.2.0, 20 November 2024
%
    arguments
        vCat double
        PP struct
        OptArgs.PointProcess char    = 'Hawkes'     % 
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
    vCat(:,2) = vCat(:,2) - PP.fM0; % subtract the reference magnitude fM0
    
    if strcmp(OptArgs.PointProcess,'Hawkes')
        gif_pp = @(t,vcat,vPar,bRateType) gif_hawkes(t,vcat,PP.fTs,PP.nJs,vPar,bRateType);
        negLLfun = @(vPar) -llfun_hawkes(vPar,vCat,PP.fTs,PP.fTe,PP.nJs,PP.nJe);
    end
    tic
    [PP.vPar, PP.fLle, exitflag, PP.mHessian] = pp_fit_mle(negLLfun,PP,'OptimMethod',OptArgs.OptimMethod);
    %disp([llfun_hawkes(PP.vPar,vCat,PP.fTs,PP.fTe,PP.nJs,PP.nJe) PP.fLle])
    PP.nPar = length(PP.vPar);
    PP.AIC = -2.0*PP.fLle + 2*PP.nPar; % standard definition
    [PP.vPar_ci, PP.vPar_err, PP.mCovMat] = pp_param_ci(PP.vPar,vCat,PP);
    toc
    disp(['Exit flag: ',num2str(exitflag)]);

    if ~isempty(OptArgs.SaveFigure)
        vCat_save = [vCat(:,1), vCat(:,2) + PP.fM0];
        save([OptArgs.SaveFigure,'_eq_cat.dat'],'vCat_save','-ascii');
        vPar_save = [PP.vPar; PP.vPar_err];
        save([OptArgs.SaveFigure,'_parameters.txt'],'vPar_save','-ascii');
    end
    
    plot_pp_rate(vCat,PP,gif_pp,...
        'PlotRate',OptArgs.PlotRate,'PlotCumulative',OptArgs.PlotCumulative,'PlotTransformed',OptArgs.PlotTransformed,'PlotCumulSum',OptArgs.PlotCumulSum,...
        'BinNumber',OptArgs.BinNumber,'DateFormat',OptArgs.DateFormat,'Title',OptArgs.Title,'SaveFigure',OptArgs.SaveFigure,'Position',OptArgs.Position);
end

