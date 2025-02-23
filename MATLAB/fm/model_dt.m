function model_dt(vCat,Model,varargin)
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 22 May 2022
%   ...
%   version 1.0.0, 22 May 2022
%
    sFitModel   = {'Exp','Gamma'}; % 'Gamma', 
    fAlpha      = 0.05;
    nBins       = 30;
    sPlotStatus = 'on';
    sBinType    = 'loglog';
    for k = 1:length(varargin)
        if strcmp('FitModel',varargin{k})
            sFitModel = varargin{k+1};
        end
        if strcmp('Alpha',varargin{k})
            fAlpha = varargin{k+1};
        end
        if strcmp('Binning',varargin{k})
            sBinType = varargin{k+1};
        end
        if strcmp('Plot',varargin{k})
            sPlotStatus = varargin{k+1};
        end
        if strcmp('NumBins',varargin{k})
            nBins = varargin{k+1};
        end
    end

    vT = vCat(:,1);
    vDt = diff(vT);
    fTmin = Model.fTmin;
    fTmax = Model.fTmax;

    if strcmp(sBinType,'linear')
        vXedges = linspace(fTmin,fTmax,nBins);
        vX = vXedges(1:end-1) + (fTmax - fTmin)/(nBins-1);
        vXmodel = linspace(fTmin,fTmax,10*nBins) + (fTmax - fTmin)/(10*nBins-1);
    elseif strcmp(sBinType,'loglog')
        vXedges = logspace(log10(fTmin),log10(fTmax),nBins); % vector of the right most coordinates of the bins
        vX = zeros(1,nBins-1);
        for k = 1:nBins-1
           vX(k) = sqrt(vXedges(k)*vXedges(k+1)); % geometric mean
        end
        vXedges_m = logspace(log10(fTmin),log10(fTmax),10*nBins); % vector of the right most coordinates of the bins
        vXmodel = zeros(1,10*nBins-1);
        for k = 1:10*nBins-1
           vXmodel(k) = sqrt(vXedges_m(k)*vXedges_m(k+1)); % geometric mean
        end
    end
    % compute the histogram
    vDt_pdf = histcounts(vDt,vXedges,'Normalization','pdf');           % histogram of the interevent times

    % limit between Tmin and Tmax
    vDt = vDt(vDt >= fTmin & vDt <= fTmax);
    nDt = length(vDt);

    if strcmp(sPlotStatus,'on')
        figure('Name','Interevent time distribution');
        %
        set(gcf,'color','w'); % background color for the figure
        % plot interevent time distribution
        if strcmp(sBinType,'linear')
            semilogy(vX,vDt_pdf,'s','MarkerEdgeColor','#EDB120','MarkerFaceColor','#EDB120','MarkerSize',8);
        elseif strcmp(sBinType,'loglog')
            loglog(vX,vDt_pdf,'s','MarkerEdgeColor','#EDB120','MarkerFaceColor','#EDB120','MarkerSize',8);
        end
        ax = gca;
        yl1 = ax.YLim(1);
        yl2 = ax.YLim(2);
        hold on
        k = 1;
        cLegend{k} = 'interevent-time distribution';
        for i = 1:length(sFitModel)
            if strcmp(sFitModel{i},'Exp')
                [muhat, muci] = expfit(vDt,fAlpha);
                plot(vXmodel,exppdf(vXmodel,muhat),'Color','#0072BD','LineWidth',1.6);
                ylim([yl1,yl2]);
                k = k + 1;
                cLegend{k} = ['exponential fit: \lambda = ',num2str(1/muhat,4),' \pm ',num2str(0.5*(muci(2)-muci(1)),2)];
                fLle = -explike(muhat,vDt);
                AIC = -2*fLle + 2;
                BIC = -2*fLle + log(nDt);
                str = ['AIC = ',num2str(AIC),', BIC = ',num2str(BIC)];
                cLegend{k} = sprintf('%s\n%s',cLegend{k},str);
            end
        end
        for i = 1:length(sFitModel)
            if strcmp(sFitModel{i},'Gamma')
                [phat, pci] = gamfit(vDt,fAlpha);
                plot(vXmodel,gampdf(vXmodel,phat(1),phat(2)),'Color','#D95319','LineWidth',1.6);
                ylim([yl1,yl2]);
                k = k + 1;
                cLegend{k} = ['gamma fit: \kappa = ',num2str(phat(1),2),' \pm ',num2str(0.5*(pci(2,1)-pci(1,1)),1),...
                              ', \theta = ',num2str(phat(2),2),' \pm ',num2str(0.5*(pci(2,2)-pci(1,2)),1)];
                fLle = -gamlike(phat,vDt);
                AIC = -2*fLle + 2*2;
                BIC = -2*fLle + 2*log(nDt);
                str = ['AIC = ',num2str(AIC),', BIC = ',num2str(BIC)];
                cLegend{k} = sprintf('%s\n%s',cLegend{k},str);
            end
        end
        legend(cLegend,'Location','southwest');
        xlabel('$\Delta t$, days','Interpreter','latex');
        ylabel('pdf','FontName','TimesNewRoman');
        title(Model.sTitle);
        hold off
        
        sName = [Model.sFileName,'_dt'];
%         vTmp = [vX', vDt_pdf'];
%         save([sName,'.dat'],'vTmp','-ascii');
        saveas(gcf,[sName,'.png']);
%        saveas(gcf,[sName,'.pdf']);
        saveas(gcf,[sName,'.eps'],'epsc');

        figure('Name','Interevent time distribution');
        %
        set(gcf,'color','w'); % background color for the figure
        co = get(0,'DefaultAxesColorOrder');
        for k = 1:length(Model.vMag_c)
            vCat_mc = vCat(vCat(:,2) >= Model.vMag_c(k),:);
            vT = vCat_mc(:,1);
            vDt = diff(vT);
            % compute the histogram
            vDt_pdf = histcounts(vDt,vXedges,'Normalization','pdf');           % histogram of the interevent times
            % plot interevent time distribution
            if strcmp(sBinType,'linear')
                semilogy(vX,vDt_pdf,'-s','MarkerEdgeColor',co(k,:),'MarkerFaceColor',co(k,:),'MarkerSize',5);
            elseif strcmp(sBinType,'loglog')
                loglog(vX,vDt_pdf,'-s','MarkerEdgeColor',co(k,:),'MarkerFaceColor',co(k,:),'MarkerSize',5);
            end
            cLegend{k} = ['m_c = ',num2str(Model.vMag_c(k))];
            hold on
        end
        legend(cLegend,'Location','southwest');
        xlabel('$\Delta t$, days','Interpreter','latex');
        ylabel('pdf','FontName','TimesNewRoman');
        title(Model.sTitle);
        hold off
        
        sName = [Model.sFileName,'_dt_mc'];
        saveas(gcf,[sName,'.png']);
%        saveas(gcf,[sName,'.pdf']);
        saveas(gcf,[sName,'.eps'],'epsc');
    end
end

