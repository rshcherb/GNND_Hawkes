function [vX, vY, vXmod, vModel, fGof] = plot_pp_gif(vCat,PP,gif_pp,OptArgs)
%
%   vCat   - the catalogue of events
%   PP     - input structure for the point process: vPar, fT0, fTs, fTe, fT1, ...
%   gif_pp - function handle for the ground intensity function of the point process 
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 31 October 2024
%   ...
%   version: 1.2.0, 3 November 2024
%
    arguments
        vCat double
        PP struct                                     % structure that defines the poit-process to fit
        gif_pp                                        % function handle for the point process ground intensity function
        OptArgs.PointProcess char = 'Hawkes'          % 
        OptArgs.DateFormat char   = 'day'             % to use for time axis as 'day' or 'date'
        OptArgs.ErrorMethod char  = 'Hessian'         % 'LogLik', 'BootStrap', 'Simul'
        OptArgs.BinNumber double  = []                % number of binns
        OptArgs.LogBinNumber double = 26              %
        OptArgs.RateType char     = 'cumulative'      % 
        OptArgs.BinType char      = 'linear'          % 'loglog'
        OptArgs.Time char         = 'original'        % 'transformed', 'cusum'
        OptArgs.ParamLegend char  = 'on'              % whether to plot the parameter legend
        OptArgs.Position double   = [400,300,800,600] % figure size and position
        OptArgs.ColorScheme char  = 'BlueScatterPurple' % 'BlueScatterOrange': blue and purple colors and scatter symbols;'OrangeStem': orange, blue colors and stem symbols
        OptArgs.sLYLabel char     = '$\lambda(t,\mathcal{H}_t)$' % 
    end
    nDataLen = length(vCat(:,1));
    if isempty(OptArgs.BinNumber)
        if nDataLen < 100
            OptArgs.BinNumber = 1000;
        else
            OptArgs.BinNumber = nDataLen;
        end
    end
    fDateS = 0;
    if strcmp(OptArgs.DateFormat,'date')
        fDateS = datenum(PP.vDateStart);
    end

    fT0 = PP.fT0;
    fTs = PP.fTs;
    fTe = PP.fTe;
    fT1 = PP.fT1;
    nPar = length(PP.vPar);
    % set the plot properties
    sposx = 0.05; % left string position
    sposy = 0.95;
    TxtFontSize = 12;
    LabelFontSize = 14;
    LgdFontSize = 14;
    MrkrFaceAlpha = 0.2;
    MrkrFaceCol = "#0072BD";
    MrkrEdgeCol = "#0072BD";
    CumulMrkrCol = "#4DBEEE";
    %CumulMrkrCol = "#77AC30";
    if strcmp(OptArgs.ColorScheme,'OrangeStem')
        MrkrFaceCol = "#D95319";
        MrkrEdgeCol = "#D95319";
        LineCol = "#0072BD";
    elseif strcmp(OptArgs.ColorScheme,'BlueScatterOrange')
        LineCol = "#D95319";
    elseif strcmp(OptArgs.ColorScheme,'BlueScatterPurple')
        %LineCol = "#000000";
        LineCol = "#7E2F8E";
    end

    vTdata = vCat(:,1);
    vMag = vCat(:,2) + PP.fM0;
    if strcmp(OptArgs.RateType,'rate')
        if strcmp(OptArgs.BinType,'linear')
            indx = (vTdata >= fTs) & (vTdata <= fTe);
            vXmod = vTdata(indx);
            %vXmod = linspace(fTs,fTe,nBinNum)'; % this gives some unusual spikes in the rate when plotted
            vTmp = [vXmod; linspace(fTs,fTe,OptArgs.BinNumber)'];
            vXmod = sort(vTmp);
            vModel = gif_pp(vXmod,vCat,PP.vPar,0);
            indx = vModel > 1e-8; % exclude zero values for the rate model 
            vModel = vModel(indx);
            vXmod = vXmod(indx);
            if strcmp(OptArgs.ColorScheme,'OrangeStem')
            elseif strcmp(OptArgs.ColorScheme,'BlueScatterOrange')
                %colororder(["#A2142F","#0072BD"]);
                colororder([LineCol,MrkrFaceCol]);
            elseif strcmp(OptArgs.ColorScheme,'BlueScatterPurple')
                colororder([LineCol,MrkrFaceCol]);
            end
            yyaxis left
            semilogy(vXmod + fDateS,vModel,'LineWidth',1.0,'Color',LineCol);
            ax = gca; ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
            ylim([0.9*min(vModel),Inf]);
            ylabel(OptArgs.sLYLabel,'Interpreter','latex','FontSize',LabelFontSize);
            if strcmp(OptArgs.ColorScheme,'OrangeStem')
                yyaxis right
                hs = stem(vTdata + fDateS,vMag,':d','BaseValue',PP.fMc,'MarkerFaceColor',MrkrFaceCol,'MarkerEdgeColor',MrkrEdgeCol);
                hbase = hs.BaseLine; hbase.Visible = 'off';
            elseif strcmp(OptArgs.ColorScheme,'BlueScatterPurple') || strcmp(OptArgs.ColorScheme,'BlueScatterOrange')
                yyaxis right
                MrkrSize = (2 + vMag - min(vMag)).^4;
                scatter(vTdata + fDateS,vMag,MrkrSize,'filled','d','MarkerFaceColor',MrkrFaceCol,'MarkerEdgeColor',MrkrEdgeCol,'MarkerFaceAlpha',MrkrFaceAlpha);
            end
            ylabel('Magnitude, $m$','Interpreter','latex','FontSize',LabelFontSize);
            ylim([min(vMag)-0.05, max(vMag)+1.0]);
            ax = gca; ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
            hold on
            plot(ax,[fTs + fDateS, fTs + fDateS],[ax.YLim(1), ax.YLim(2)],'k:',...
                    [fTe + fDateS, fTe + fDateS],[ax.YLim(1), ax.YLim(2)],'k:');
            hold off
            axis tight
            text(fTs + fDateS,0.99*ax.YLim(2),'$T_s$','Interpreter','latex','FontSize',TxtFontSize);
            text(fTe + fDateS,0.99*ax.YLim(2),'$T_e$','Interpreter','latex','FontSize',TxtFontSize);
            vX = vCat(:,1);
            vY = vMag;
            if strcmp(OptArgs.DateFormat,'day')
                xlabel('$t$','Interpreter','latex','FontSize',LabelFontSize);
            elseif strcmp(OptArgs.DateFormat,'date')
                xlabel('Date (mm/yyyy)','Interpreter','latex','FontSize',LabelFontSize);
                %datetick('x','mm/yy','keepticks');
                datetick('x','mm/yyyy','keeplimits');
            end
            % this does not work properly with MarkerFaceAlpha propety
            % get(gca,'SortMethod');
            % set(gca,'SortMethod','depth'); % makes the left axis on top of the right axis: https://stackoverflow.com/questions/37709988/plot-line-chart-on-the-left-axis-uper-the-bar-chart-on-the-right-axis
            % lp.ZData = ones(size(lp.XData));
            % sp.ZData = ones(size(sp.XData));
            %
            % another possibility: https://www.mathworks.com/matlabcentral/answers/515212-when-using-yyaxis-is-there-a-way-to-plot-the-right-yaxis-data-behind-the-left-yaxis-data
        elseif strcmp(OptArgs.BinType,'loglog')
            [vData, vXmod, vBin] = event_rate(vTdata,'Binning','loglog',OptArgs.LogBinNumber);
            vModel = gif_pp(vXmod,vCat,PP.vPar,0);
            pl1 = loglog(vXmod,vModel,'-','LineWidth',1.5,'Color','#7E2F8E');
            xlim([fTs,fTe]);
            grid on
            hold on
            pl2 = loglog(vXmod,vData,'s','MarkerSize',8,'MarkerFaceColor','#4DBEEE','MarkerEdgeColor','#4DBEEE');
            ylabel('$\lambda(t)$','Interpreter','latex','FontSize',LabelFontSize);
            vX = vXmod;
            vY = vData;
            sposy = 0.35;
            ax = gca; ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
            plot(ax,[fTs, fTs],[ax.YLim(1), ax.YLim(2)],'k:',[fTe, fTe],[ax.YLim(1), ax.YLim(2)],'k:');
            hold off
            text(fTs,0.99*ax.YLim(2),'$T_s$','Interpreter','latex','FontSize',TxtFontSize);
            text(fTe,0.99*ax.YLim(2),'$T_e$','Interpreter','latex','FontSize',TxtFontSize);
            xlabel('$t$','Interpreter','latex','FontSize',LabelFontSize);
            legend([pl1, pl2],'model','observed','Interpreter','latex','FontSize','FontSize',LgdFontSize);
        end
    elseif strcmp(OptArgs.RateType,'cumulative')
        vCumData = (1:nDataLen)';
        fCumData0 = find(vTdata(:,1) < fTs,1,'last'); % find the last index for which vTdata < fTs
        if isempty(fCumData0)
            fCumData0 = 0;
        end
        
        % original time: we define nBinNum number of points to compute the rate
        if strcmp(OptArgs.Time,'original')
            vT = linspace(fTs,fTe,OptArgs.BinNumber)';
        elseif strcmp(OptArgs.Time,'transformed') || strcmp(OptArgs.Time,'cusum')
            indx = (vTdata >= fTs) & (vTdata <= fTe);
            vT = vTdata(indx);
        end
        vCumModel = fCumData0 + gif_pp(vT,vCat,PP.vPar,1);
        
        if strcmp(OptArgs.Time,'original')
            vX = vTdata + fDateS;  % x data
            vY = vCumData;         % y data
            vXmod = vT + fDateS;
            vModel = vCumModel;
        elseif strcmp(OptArgs.Time,'transformed')
            vX = vCumModel;
            vY = vCumData(PP.nJs:PP.nJe);
            vXmod = vCumModel;
            vModel = vCumModel;
        elseif strcmp(OptArgs.Time,'cusum')
            vX = vCumModel;
            vY = vCumData(PP.nJs:PP.nJe) - vCumModel;
            vXmod = vCumModel;
            vModel = zeros(length(vT),1);
        end
        
        if strcmp(OptArgs.Time,'original') || strcmp(OptArgs.Time,'transformed')
            plot(vX,vY,'.','Color',CumulMrkrCol,'MarkerSize',12); %
        elseif strcmp(OptArgs.Time,'cusum')
            plot(vX,vY,'.:','Color',CumulMrkrCol); % 
        end
        ax = gca; ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
        axis tight
        grid on
        hold on
        plot(vXmod,vModel,'-','LineWidth',1.6,'Color',LineCol); % 
        if strcmp(OptArgs.Time,'original')
            hold on
            plot([fTs + fDateS, fTs + fDateS],[ax.YLim(1), ax.YLim(2)],'k:',...
                 [fTe + fDateS, fTe + fDateS],[ax.YLim(1), ax.YLim(2)],'k:');
            text(fTs + fDateS,0.99*ax.YLim(2),'$T_s$','Interpreter','latex','FontSize',TxtFontSize);
            text(fTe + fDateS,0.99*ax.YLim(2),'$T_e$','Interpreter','latex','FontSize',TxtFontSize);
            if strcmp(OptArgs.DateFormat,'day')
                xlabel('$t$','Interpreter','latex','FontSize',LabelFontSize);
            elseif strcmp(OptArgs.DateFormat,'date')
                xlabel('date (mm/yyyy)','Interpreter','latex','FontSize',LabelFontSize);
                %datetick('x','mm/yy','keepticks');
                datetick('x','mm/yyyy','keeplimits');
            end
            ylabel('Cumulative number, $N$','Interpreter','latex','FontSize',LabelFontSize);
            legend('observed','model','Interpreter','latex','Location','southeast','FontSize',LgdFontSize);
        elseif strcmp(OptArgs.Time,'transformed')
            xlabel('Transformed time','Interpreter','latex','FontSize',LabelFontSize);
            ylabel('Event number','Interpreter','latex','FontSize',LabelFontSize);
            legend('observed','model','Interpreter','latex','Location','southeast','FontSize',LgdFontSize);
        elseif strcmp(OptArgs.Time,'cusum')
            xlabel('Transformed time','Interpreter','latex','FontSize',LabelFontSize);
            ylabel('cusum','Interpreter','latex','FontSize',LabelFontSize);
        end
        hold off
    end
    fGof = pp_gof(vCat,PP,gif_pp);
    if ~(strcmp(OptArgs.RateType,'rate') && strcmp(OptArgs.BinType,'linear'))
        % add parameter legend to cumulative plots
        if strcmp(OptArgs.ParamLegend,'on')
            dy = 0.05;
            for n = 1:nPar
                %str = sprintf('%s = %.4g \xB1 %.4g',vParName{n},vPar(n),vParErr(n));
                str = ['$',PP.ParName{n},' = ',num2str(PP.vPar(n),'%.3g'),' \pm ',num2str(PP.vPar_err(n),'%.3g'),'$'];
                text(sposx,sposy-(n-1)*dy,str,'Units','normal','Interpreter','latex','FontSize',TxtFontSize);        
            end
            str = sprintf('LL = %.3f; AIC = %.4f',PP.fLle,PP.AIC);
            text(sposx,sposy-nPar*dy,str,'Units','normal','Interpreter','latex','FontSize',TxtFontSize);
            str = sprintf('gof = %.4f',fGof);
            text(sposx,sposy-(nPar+1)*dy,str,'Units','normal','Interpreter','latex','FontSize',TxtFontSize);
        end
    end
end

