function plot_pp_rate(vCat,PP,gif_pp,OptArgs)
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 31 October 2024
%   ...
%   version: 1.1.0, 20 November 2024
%
    arguments
        vCat double
        PP struct                                   % structure that defines the poit-process to fit
        gif_pp                                      % function handle for the point process ground intensity function
        OptArgs.PointProcess char    = 'Hawkes'     % 'ETAS'
        OptArgs.DateFormat char      = 'day'        % to use for time axis as 'day' or 'date'
        OptArgs.ErrorMethod char     = 'Hessian'    % 'LogLik', 'BootStrap', 'Simul'
        OptArgs.BinNumber double     = []           % number of binns
        OptArgs.PlotRate char        = 'on'         % 
        OptArgs.PlotCumulative char  = 'on'         % 
        OptArgs.PlotTransformed char = 'on'         % 
        OptArgs.PlotCumulSum char    = 'on'         % 
        OptArgs.Title char           = []           %  
        OptArgs.SaveFigure char      = []           % 
        OptArgs.Position double = [400,300,800,600] % figure size and position
    end

    if OptArgs.PlotRate
        figure('Name','Hawkes: rate','Position',OptArgs.Position);
        set(gcf,'color','w'); % background color for the figure
        [vX, vY, vXmod, vModel, fGof] = plot_pp_gif(vCat,PP,gif_pp,'RateType','rate','BinNumber',OptArgs.BinNumber,'DateFormat',OptArgs.DateFormat);
        title(OptArgs.Title);
        if ~isempty(OptArgs.SaveFigure)
            sName = [OptArgs.SaveFigure,'_rate'];
            save_cf(gcf,sName,'fig','png','pdf');
            mTmp = [vXmod, vModel];
            save([sName,'.dat'],'mTmp','-ascii');
        end
    end
    
    if OptArgs.PlotCumulative
        figure('Name','Hawkes: cumulative','Position',OptArgs.Position);
        set(gcf,'color','w'); % background color for the figure
        [vX, vY, vXmod, vModel, fGof] = plot_pp_gif(vCat,PP,gif_pp,'BinNumber',OptArgs.BinNumber,'DateFormat',OptArgs.DateFormat);
        title(OptArgs.Title);
        if ~isempty(OptArgs.SaveFigure)
            sName = [OptArgs.SaveFigure,'_cumulative'];
            save_cf(gcf,sName,'fig','png','pdf');
            mTmp = [vXmod, vModel];
            save([sName,'.dat'],'mTmp','-ascii');
            mTmp = [vX, vY];
            save([OptArgs.SaveFigure,'_eq_cumulative.dat'],'mTmp','-ascii');
        end
    end
    
    if OptArgs.PlotTransformed
        figure('Name','Hawkes: cumulative - transformed time','Position',OptArgs.Position);
        set(gcf,'color','w'); % background color for the figure
        [vX, vY, vXmod, vModel, fGof] = plot_pp_gif(vCat,PP,gif_pp,'Time','transformed','DateFormat',OptArgs.DateFormat);
        title(OptArgs.Title);
        if ~isempty(OptArgs.SaveFigure)
            sName = [OptArgs.SaveFigure,'_transformed'];
            save_cf(gcf,sName,'fig','png','pdf');
            mTmp = [vX, vY];
            save([sName,'.dat'],'mTmp','-ascii');
        end
    end
    
    if OptArgs.PlotCumulSum
        figure('Name','Hawkes: cumulative - cusum time','Position',OptArgs.Position);
        set(gcf,'color','w'); % background color for the figure
        [vX, vY, vXmod, vModel, fGof] = plot_pp_gif(vCat,PP,gif_pp,'Time','cusum','DateFormat',OptArgs.DateFormat);
        title(OptArgs.Title);
        if ~isempty(OptArgs.SaveFigure)
            sName = [OptArgs.SaveFigure,'_cusum'];
            save_cf(gcf,sName,'fig','png','pdf');
            mTmp = [vX, vY];
            save([sName,'.dat'],'mTmp','-ascii');
        end
    end
end

