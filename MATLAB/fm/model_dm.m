function model_dm(vMag,Model,varargin)
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 23 May 2022
%   ...
%   version 1.1.0, 1 December 2022
%
    fAlpha      = 0.05;
    nBins       = 101;
    sPlotStatus = 'on';
    for k = 1:length(varargin)
        if strcmp('Alpha',varargin{k})
            fAlpha = varargin{k+1};
        end
        if strcmp('Plot',varargin{k})
            sPlotStatus = varargin{k+1};
        end
        if strcmp('NumBins',varargin{k})
            nBins = varargin{k+1};
        end
    end
%     nMag = length(vMag);
%     vMag = vMag(randperm(nMag));
    nSampl = Model.nSampl;
    %[vX, vDDm_cdf] = comp_dcdf_dm(vMag,Model.fDMmin,Model.fDMmax,nBins,Model.fMc);
    vDDm_cdf = zeros(1,nBins-1);
    for n = 1:nSampl
        [vX, vDDm_cdf_k] = comp_dcdf_dm(vMag,Model.fDMmin,Model.fDMmax,nBins,Model.fMc);
        vDDm_cdf = vDDm_cdf + vDDm_cdf_k;
    end
    vDDm_cdf = vDDm_cdf/double(nSampl);
    vStdErr = comp_errors(vMag,Model.fDMmin,Model.fDMmax,Model.nSampl,nBins,Model.fSigma,Model.fMc);

    if strcmp(sPlotStatus,'on')
        figure('Name','Distribution of magnitude differences','Position',[550 400 800 600]);
        %
        set(gcf,'color','w'); % background color for the figure
        % plot the distribution of magnitude differences
        plot(vX,vDDm_cdf,'-o','Color','#EDB120','MarkerEdgeColor','#EDB120','MarkerFaceColor','#EDB120','MarkerSize',4);
        ax = gca; ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
        hold on
        errorbar(vX,vDDm_cdf,vStdErr,'k','LineStyle','none',...
            'Marker','o','MarkerEdgeColor','#EDB120','MarkerFaceColor','#EDB120','MarkerSize',4);

        legend({'$\Delta P(dm)$',[num2str(Model.fSigma,1),'$\sigma$ error']},'Location','northwest','Interpreter','latex');
        xlabel('$dm$','Interpreter','latex');
        ylabel('$\Delta P(dm)$','Interpreter','latex');
        title(Model.sTitle);
        hold off
        
        sName = [Model.sFileName,'_dm'];
        vTmp = [vX', vDDm_cdf',vStdErr'];
        save([sName,'.dat'],'vTmp','-ascii');
        saveas(gcf,[sName,'.png']);
        saveas(gcf,[sName,'.pdf']);
        saveas(gcf,[sName,'.fig']);
        saveas(gcf,[sName,'.eps'],'epsc');

        figure('Name','Distribution of magnitude differences with changing m_c','Position',[550 400 800 600]);
        %
        set(gcf,'color','w'); % background color for the figure
        co = get(0,'DefaultAxesColorOrder');
        for k = 1:length(Model.vMag_c)
            %[vX, vDDm_cdf] = comp_dcdf_dm(vMag,Model.fDMmin,Model.fDMmax,nBins,Model.vMag_c(k));
            vDDm_cdf = zeros(1,nBins-1);
            for n = 1:nSampl
                [vX, vDDm_cdf_k] = comp_dcdf_dm(vMag,Model.fDMmin,Model.fDMmax,nBins,Model.vMag_c(k));
                vDDm_cdf = vDDm_cdf + vDDm_cdf_k;
            end
            vDDm_cdf = vDDm_cdf/double(nSampl);
            
            % plot the distribution of magnitude differences
            %plot(vX,vDDm_cdf,'-o','Color','#EDB120','MarkerEdgeColor','#EDB120','MarkerFaceColor','#EDB120','MarkerSize',4);
            plot(vX,vDDm_cdf,'-o','Color',co(k,:),'MarkerEdgeColor',co(k,:),'MarkerFaceColor',co(k,:),'MarkerSize',4);
            %pl = plot(vX,vDDm_cdf,'-o','MarkerFaceColor','auto','MarkerSize',4);
            ax = gca; ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
            cLegend{k} = ['m_c = ',num2str(Model.vMag_c(k))];
            hold on
        end
        legend(cLegend,'Location','northwest');
        xlabel('$dm$','Interpreter','latex');
        ylabel('$\Delta P(dm)$','Interpreter','latex');
        title(Model.sTitle);
        hold off
        
        sName = [Model.sFileName,'_dm_mc'];
        saveas(gcf,[sName,'.png']);
        saveas(gcf,[sName,'.fig']);
        saveas(gcf,[sName,'.pdf']);
        saveas(gcf,[sName,'.eps'],'epsc');
    end
end

%
function [vX, vDDm_cdf] = comp_dcdf_dm(vMag,fDMmin,fDMmax,nBins,fMc)
%
    vMag       = vMag(vMag >= fMc);
    vXedges    = linspace(fDMmin,fDMmax,nBins);
    vX         = vXedges(1:end-1) + (fDMmax - fDMmin)/(nBins-1);
    nMag       = length(vMag);
    %vMag = vMag(randperm(nMag));                         % reshuffled magnitudes
    vDm        = diff(vMag);                              % magnitude difference between the original magnitudes
    vMag_shufl = vMag(randperm(nMag));                    % reshuffled magnitudes
    vDm_shufl  = diff(vMag_shufl);
    nDMag      = length(vDm);

    % compute the histogram
    vDm_hist       = histcounts(vDm,vXedges);             % histogram of the magnitude differences
    vDm_shufl_hist = histcounts(vDm_shufl,vXedges);       % histogram of the magnitude differences
    vDm_cdf        = cumsum(vDm_hist)/double(nDMag);
    vDm_shufl_cdf  = cumsum(vDm_shufl_hist)/double(nDMag);
    vDDm_cdf       = vDm_cdf - vDm_shufl_cdf;
end

%
function vStdErr = comp_errors(vMag,fDMmin,fDMmax,nSampl,nBins,fSigma,fMc)
%

    vDDm_cdf_mean = zeros(1,nBins-1);
    vDDm_cdf_var = zeros(1,nBins-1);
    for n = 1:nSampl
        [vX, vDDm_cdf] = comp_dcdf_dm(vMag,fDMmin,fDMmax,nBins,fMc);
        % compute the difference betweeen cumulative distributions
        vDDm_cdf_mean = vDDm_cdf_mean + vDDm_cdf;
        vDDm_cdf_var = vDDm_cdf_var + vDDm_cdf.^2;
    end

    vDDm_cdf_mean = vDDm_cdf_mean/double(nSampl);
    vDDm_cdf_var  = vDDm_cdf_var/double(nSampl);
    dtmp          = vDDm_cdf_var - vDDm_cdf_mean.^2;
    vStdErr       = fSigma*sqrt(double(nSampl)/(double(nSampl)-1.0).*dtmp);

%     vDDm_cdf_var = zeros(1,nBins-1);
%     for n = 1:nSampl
%         [vX, vDDm_cdf] = comp_dcdf_dm(vMag,fDMmin,fDMmax,nBins,fMc);
%         % compute the difference betweeen cumulative distributions
%         %vDDm_cdf_mean = vDDm_cdf_mean + vDDm_cdf;
%         vDDm_cdf_var = vDDm_cdf_var + (vDDm_cdf - vDDm_cdf_mean).^2;
%     end
%     vDDm_cdf_var  = vDDm_cdf_var/(double(nSampl) - 1.0);
%     vStdErr = fSigma*sqrt(vDDm_cdf_var);
end

