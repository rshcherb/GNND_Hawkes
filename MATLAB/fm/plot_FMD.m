function [vX, vHist_emp, vPdf_emp, vCdf_emp] = plot_FMD(vMag,fMc,fDm,OptArgs)
%
%   Function to plot the frequency-momnet statistics
%
%   vMag  - list of magnitudes to bin and construct the distribution
%   fMc   - lower magnitude cutoff
%   fDm   - the binning of maghnitudes
%
%   Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 17 June 2022
%   ...
%   version: 3.1.0, 11 January 2025
%
    arguments
        vMag double                    % earthquake magnitudes
        fMc double                     % lower magnitude cutoff
        fDm double                     % magnitude binning
        OptArgs.Units char = 'moment'  % in which units to plot: 'moment' or 'magnitude'
        OptArgs.Pareto struct = []     % structure for Pareto model parameters
        OptArgs.TrnPareto struct = []  % structure for truncated Pareto model parameters
        OptArgs.TapPareto struct = []  % structure for tapered Pareto model parameters
        OptArgs.TapTrnPareto struct = []  % structure for tapered Pareto + truncated Pareto model parameters
        OptArgs.TapTapPareto struct = []  % structure for tapered Pareto + tapered Pareto model parameters
        OptArgs.TapParetoPareto struct = []  % structure for tapered Pareto + Pareto model parameters
        OptArgs.ParetoPareto struct = []  % structure for Pareto + Pareto model parameters
        OptArgs.Xmin double = []       % lower cutoff
        OptArgs.Xmax double = []       % upper cutoff
        OptArgs.Alpha double           % confidence level
    end
    if isempty(OptArgs.Xmin)
        %OptArgs.Xmin = floor(min(vMag)/fDm)*fDm;
        OptArgs.Xmin = fMc;
    end
    if isempty(OptArgs.Xmax)
        OptArgs.Xmax = ceil(max(vMag)/fDm)*fDm; % fix the upper limit
    end
    
    M_min = mag2moment(fMc);
    if strcmp(OptArgs.Units,'moment')
        nBins = round((OptArgs.Xmax - OptArgs.Xmin)/fDm/2);
        vMom  = mag2moment(vMag);   % convert to moment
        fMmin = mag2moment(OptArgs.Xmin);  % convert to moment
        fMmax = mag2moment(OptArgs.Xmax);  % convert to moment
        a1 = log10(fMmin);
        a2 = log10(fMmax);
        vXedges = logspace(a1,a2,nBins+1);
        vX = zeros(1,nBins);
        for n = 1:nBins
            vX(n) = geomean([vXedges(n),vXedges(n+1)]); % take the geometric mean
        end
        vPdf_mom = histcounts(vMom,vXedges,'Normalization','pdf');
        vHist_mom = vPdf_mom;
        vCdf_mom = [];
        vX_mom = logspace(a1,a2,300);
        vHist_emp = vHist_mom; vPdf_emp = vPdf_mom; vCdf_emp = vCdf_mom;
    end
    if strcmp(OptArgs.Units,'magnitude')
        vX = OptArgs.Xmin:fDm:OptArgs.Xmax;   % centres of the bins, total N, between Xmin and Xmax 
        %                                                            vX(1)        vX(k)
        vXedges = (OptArgs.Xmin-0.5*fDm):fDm:(OptArgs.Xmax+0.5*fDm); % bin edges: |_____|_..._|_____|_..._|_____|, total N+1
        %                                                           1     2     k    k+1
        % compute the histogram
        vHist_mag = histcounts(vMag,vXedges);           % histogram of the earthquake magnitudes
        nBins = length(vHist_mag);                       % number of bins, N
        vPdf_mag = histcounts(vMag,vXedges,'Normalization','pdf');
        % complimentary cumulative distribution
        vCdf_mag = [];
        %indx = vHist_mag > 0; % use only nonzero bins
        %vHist_mag = vHist_mag(indx);
        %vPdf_mag = vPdf_mag(indx);
        %vX = vX(indx);
        vXmag = linspace(OptArgs.Xmin,OptArgs.Xmax,300);
        vX_mom = mag2moment(vXmag);
        vHist_emp = vHist_mag; vPdf_emp = vPdf_mag; vCdf_emp = vCdf_mag;
    end

    lwidth       = 1.6;       % line width
    mrkrcol      = '#4DBEEE'; % marker color
    mrkrsize     = 8;         % marker size
    colpar       = '#0072BD'; % Pareto
    coltrnpar    = '#EDB120'; % truncated Pareto
    coltappar    = '#77AC30'; % tapered Pareto
    coltaptrnpar = '#A2142F'; % tapered Pareto + truncated Pareto
    coltaptappar = '#D95319'; % tapered Pareto + tapered Pareto
    coltapparpar = '#7E2F8E'; % tapered Pareto + Pareto
    colparpar    = '#0000FF'; % Pareto + Pareto
    set(gcf,'color','w'); % background color for the figure
    if strcmp(OptArgs.Units,'moment')
        tl = tiledlayout(1,1);
        ax = axes(tl);
        loglog(ax,vX,vPdf_mom,'s','MarkerEdgeColor','none','MarkerFaceColor',mrkrcol,'MarkerSize',mrkrsize);
        xlabel(ax,'Moment, $M$','Interpreter','latex');
        ax.Box = 'off';
        ymin = 0.8*min(vPdf_mom(vPdf_mom > 0.0));
        ymax = 1.2*max(vPdf_mom);
        %ylim([ymin, ymax]);
    end
    if strcmp(OptArgs.Units,'magnitude')
        semilogy(vX,vPdf_mag,'s','MarkerEdgeColor','none','MarkerFaceColor',mrkrcol,'MarkerSize',mrkrsize);
        ax = gca;
        xlabel('Magnitude, $m$','Interpreter','latex');
        set(gca,'XMinorTick','on');
        ymin = 0.9*min(vPdf_mag(vPdf_mag > 0.0));
        ymax = 1.1*max(vPdf_mag);
    end
    axis(ax,'manual');
    %axis tight
    grid on
    set(ax,'Layer','top');
    ylabel(ax,'pdf','Interpreter','latex');
    %ylim('manual');
    %ylim([ymin, ymax]);
        
    il = 1;
    cLegend{il} = 'observed data';
    hold on
    if ~isempty(OptArgs.Pareto)
        vPdf = pareto_pdf(vX_mom,OptArgs.Pareto.par(1),M_min); % Kagan (GJI 2002 p.521)
        if strcmp(OptArgs.Units,'moment')
            loglog(ax,vX_mom,vPdf,'--','Color',colpar,'LineWidth',lwidth);
            AIC = OptArgs.Pareto.AICmom;
        end
        if strcmp(OptArgs.Units,'magnitude')
            vPdf_mag = pdf_moment2mag(vXmag,vPdf);
            semilogy(vXmag,vPdf_mag,'--','Color',colpar,'LineWidth',lwidth);
            AIC = OptArgs.Pareto.AICmag;
        end
        il = il + 1;
        cLegend{il} = ['Pareto: $\beta = ',num2str(OptArgs.Pareto.par(1),'%.2f'),'\pm',num2str(OptArgs.Pareto.err(1),'%.2f'),...
                       '$, gof = ',num2str(OptArgs.Pareto.gof,'%.1f'),'\%, AIC = ',num2str(AIC,'%.3f')];
    end
    if ~isempty(OptArgs.TrnPareto)
        gamma = OptArgs.TrnPareto.par(1);
        M_max = OptArgs.TrnPareto.par(2);
        vPdf = trn_pareto_pdf(vX_mom,gamma,M_min,M_max); % Kagan (GJI 2002 p.523)
        strp1 = ['trn.Pareto: $\gamma = ',num2str(OptArgs.TrnPareto.par(1),'%.2f'),'\pm',num2str(OptArgs.TrnPareto.err(1),'%.2f'),'$, '];
        if strcmp(OptArgs.Units,'moment')
            loglog(ax,vX_mom,vPdf,'--','Color',coltrnpar,'LineWidth',lwidth);
            strp2 = ['$M_{max} = $',sprintf('%.2e',OptArgs.TrnPareto.par(2)),'$\,\pm\,$',sprintf('%.2e',OptArgs.TrnPareto.err(2))];
            AIC = OptArgs.TrnPareto.AICmom;
        end
        if strcmp(OptArgs.Units,'magnitude')
            vPdf_mag = pdf_moment2mag(vXmag,vPdf);
            semilogy(vXmag,vPdf_mag,'-','Color',coltrnpar,'LineWidth',lwidth);
            strp2 = ['$m_{max} = $',num2str(moment2mag(M_max),'%.2f'),'$\,\pm\,$',sprintf('%.2f',OptArgs.TrnPareto.magerr(1))];
            AIC = OptArgs.TrnPareto.AICmag;
        end
        il = il + 1;
        str2 = ['gof = ',num2str(OptArgs.TrnPareto.gof,'%.1f'),'\%, AIC = ',num2str(AIC,'%.3f')];
        cLegend{il} = sprintf('%s%s\n%s',strp1,strp2,str2);
    end
    if ~isempty(OptArgs.TapPareto)
        alpha = OptArgs.TapPareto.par(1);
        M_cm  = OptArgs.TapPareto.par(2);
        vPdf = tap_pareto_pdf(vX_mom,alpha,M_min,M_cm); % Kagan (GJI 2002 p.523)
        strp1 = ['tap.Pareto: $\alpha = ',num2str(OptArgs.TapPareto.par(1),'%.2f'),'\pm',num2str(OptArgs.TapPareto.err(1),'%.2f'),'$, '];
        if strcmp(OptArgs.Units,'moment')
            loglog(ax,vX_mom,vPdf,'-','Color',coltappar,'LineWidth',lwidth);
            strp2 = ['$M_{cm} = $ ',sprintf('%.2e',OptArgs.TapPareto.par(2)),'$\,\pm\,$',sprintf('%.2e',OptArgs.TapPareto.err(2))];
            AIC = OptArgs.TapPareto.AICmom;
        end
        if strcmp(OptArgs.Units,'magnitude')
            vPdf_mag = pdf_moment2mag(vXmag,vPdf);
            semilogy(vXmag,vPdf_mag,'-','Color',coltappar,'LineWidth',lwidth);
            strp2 = ['$m_{cm} = $ ',num2str(moment2mag(M_cm),'%.2f'),'$\,\pm\,$',sprintf('%.2f',OptArgs.TapPareto.magerr(1))];
            AIC = OptArgs.TapPareto.AICmag;
        end
        il = il + 1;
        str2 = ['gof = ',num2str(OptArgs.TapPareto.gof,'%.1f'),'\%, AIC = ',num2str(AIC,'%.3f')];
        cLegend{il} = sprintf('%s%s\n%s',strp1,strp2,str2);
    end
    if ~isempty(OptArgs.TapTrnPareto)
        alpha = OptArgs.TapTrnPareto.par(1);
        M_cm  = mag2moment(OptArgs.TapTrnPareto.par(2));
        gamma = OptArgs.TapTrnPareto.par(3);
        M_max = mag2moment(OptArgs.TapTrnPareto.par(4));
        w1 = OptArgs.TapTrnPareto.par(5);
        w2 = OptArgs.TapTrnPareto.oar(6);
        vTap = tap_pareto_pdf(vX_mom,alpha,M_min,M_cm); % Kagan (GJI 2002 p.523)
        vTrn = trn_pareto_pdf(vX_mom,gamma,M_min,M_max);
        vPdf = w1*vTap + w2*vTrn;
        strp1 = ['tap.+trn.Pareto: $\alpha = ',num2str(OptArgs.TapTrnPareto.par(1),'%.2f'),'\pm',num2str(OptArgs.TapTrnPareto.err(1),'%.2f'),'$, '];
        strp3 = ['$\gamma = ',num2str(OptArgs.TapTrnPareto.par(3),'%.2f'),'\pm',num2str(OptArgs.TapTrnPareto.err(3),'%.2f'),'$, '];
        if strcmp(OptArgs.Units,'moment')
            loglog(ax,vX_mom,vPdf,'-.','Color',coltaptrnpar,'LineWidth',lwidth);
            AIC = OptArgs.TapTrnPareto.AICmom;
        end
        if strcmp(OptArgs.Units,'magnitude')
            AIC = OptArgs.TapTrnPareto.AICmag;
        end
        il = il + 1;
        str2 = ['gof = ',num2str(OptArgs.TapTrnPareto.gof,'%.1f'),'\%, AIC = ',num2str(AIC,'%.3f')];
        cLegend{il} = sprintf('%s%s\n%s',strp1,strp3,str2);
    end
    if ~isempty(OptArgs.TapTapPareto)
        alpha1 = OptArgs.TapTapPareto.par(1);
        M_cm1  = OptArgs.TapTapPareto.par(2);
        alpha2 = OptArgs.TapTapPareto.par(3);
        M_cm2  = OptArgs.TapTapPareto.par(4);
        w      = OptArgs.TapTapPareto.par(5);
        vPdf = taptap_pareto_pdf(vX_mom,M_min,alpha1,M_cm1,alpha2,M_cm2,w);
        strp1 = ['tap.Pareto+tap.Pareto: $\alpha_1 = ',num2str(OptArgs.TapTapPareto.par(1),'%.2f'),'\pm',num2str(OptArgs.TapTapPareto.err(1),'%.2f'),'$, '];
        strp3 = ['$\alpha_2 = ',num2str(OptArgs.TapTapPareto.par(3),'%.2f'),'\pm',num2str(OptArgs.TapTapPareto.err(3),'%.2f'),'$, '];
        strp5 = ['$w = ',num2str(OptArgs.TapTapPareto.par(5),'%.2f'),'\pm',num2str(OptArgs.TapTapPareto.err(5),'%.2f'),'$, '];
        if strcmp(OptArgs.Units,'moment')
            loglog(ax,vX_mom,vPdf,'-','Color',coltaptappar,'LineWidth',lwidth);
            strp2 = ['$M_{cm}^1 = $ ',sprintf('%.2e',OptArgs.TapTapPareto.par(2)),'$\,\pm\,$',sprintf('%.2e',OptArgs.TapTapPareto.err(2)),', '];
            strp4 = ['$M_{cm}^2 = $ ',sprintf('%.2e',OptArgs.TapTapPareto.par(4)),'$\,\pm\,$',sprintf('%.2e',OptArgs.TapTapPareto.err(4)),', '];
            AIC = num2str(OptArgs.TapTapPareto.AICmom);
        end
        if strcmp(OptArgs.Units,'magnitude')
            vPdf_mag = pdf_moment2mag(vXmag,vPdf);
            semilogy(vXmag,vPdf_mag,'-','Color',coltaptappar,'LineWidth',lwidth);
            strp2 = ['$m_{cm}^1 = $ ',num2str(moment2mag(M_cm1),'%.2f'),'$\,\pm\,$',sprintf('%.2f',OptArgs.TapTapPareto.magerr(1))];
            strp4 = ['$m_{cm}^2 = $ ',num2str(moment2mag(M_cm2),'%.2f'),'$\,\pm\,$',sprintf('%.2f',OptArgs.TapTapPareto.magerr(2))];
            AIC = num2str(OptArgs.TapTapPareto.AICmag);
        end
        il = il + 1;
        str = ['gof = ',num2str(OptArgs.TapTapPareto.gof,'%.1f'),'\%, AIC = ',num2str(AIC,'%.3f')];
        cLegend{il} = sprintf('%s%s\n%s%s\n%s%s',strp1,strp2,strp3,strp4,strp5,str);
    end
    if ~isempty(OptArgs.TapParetoPareto)
        alpha = OptArgs.TapParetoPareto.par(1);
        M_cm  = OptArgs.TapParetoPareto.par(2);
        beta  = OptArgs.TapParetoPareto.par(3);
        w     = OptArgs.TapParetoPareto.par(4);
        vPdf = tap_pareto_pareto_pdf(vX_mom,M_min,alpha,M_cm,beta,w);
        strp1 = ['tap.Pareto+Pareto: $\alpha = ',num2str(OptArgs.TapParetoPareto.par(1),'%.2f'),'\pm',num2str(OptArgs.TapParetoPareto.err(1),'%.2f'),'$, '];
        strp3 = ['$\beta = ',num2str(OptArgs.TapParetoPareto.par(3),'%.2f'),'\pm',num2str(OptArgs.TapParetoPareto.err(3),'%.2f'),'$, '];
        strp4 = ['$w = ',num2str(OptArgs.TapParetoPareto.par(4),'%.2f'),'\pm',num2str(OptArgs.TapParetoPareto.err(4),'%.2f'),'$, '];
        if strcmp(OptArgs.Units,'moment')
            loglog(ax,vX_mom,vPdf,'-.','Color',coltapparpar,'LineWidth',lwidth);
            strp2 = ['$M_{cm} = $ ',sprintf('%.2e',OptArgs.TapParetoPareto.par(2)),'$\,\pm\,$',sprintf('%.2e',OptArgs.TapParetoPareto.err(2)),', '];
            AIC = OptArgs.TapParetoPareto.AICmom;
        end
        if strcmp(OptArgs.Units,'magnitude')
            vPdf_mag = pdf_moment2mag(vXmag,vPdf);
            semilogy(vXmag,vPdf_mag,'-.','Color',coltapparpar,'LineWidth',lwidth);
            strp2 = ['$m_{cm} = $ ',num2str(moment2mag(M_cm),'%.2f'),'$\,\pm\,$',sprintf('%.2f',OptArgs.TapParetoPareto.magerr(1))];
            AIC = OptArgs.TapParetoPareto.AICmag;
        end
        il = il + 1;
        str = ['gof = ',num2str(OptArgs.TapParetoPareto.gof,'%.1f'),'\%, AIC = ',num2str(AIC,'%.3f')];
        cLegend{il} = sprintf('%s%s\n%s%s%s',strp1,strp2,strp3,strp4,str);
    end
    if ~isempty(OptArgs.ParetoPareto)
        gamma = OptArgs.ParetoPareto.par(1);
        M_cm  = OptArgs.ParetoPareto.par(2);
        beta  = OptArgs.ParetoPareto.par(3);
        vPdf = pareto_pareto_pdf(vX_mom,M_min,gamma,M_cm,beta);
        if strcmp(OptArgs.Units,'moment')
            loglog(ax,vX_mom,vPdf,'-','Color',colparpar,'LineWidth',lwidth);
            AIC = OptArgs.ParetoPareto.AICmom;
        end
        if strcmp(OptArgs.Units,'magnitude')
            AIC = OptArgs.ParetoPareto.AICmag;
        end
        il = il + 1;
        cLegend{il} = ['Pareto+Pareto: gof = ',num2str(OptArgs.ParetoPareto.gof,'%.1f'),'\%, AIC = ',num2str(AIC,'%.3f')];
    end
    ax.YLim = [ymin, ymax];
    %disp(ax.YLim)
    legend(gca,cLegend,'Interpreter','latex','Location','southwest');
    hold off
    if strcmp(OptArgs.Units,'moment')
        %ax2 = axes(gcf,'Position',ax1.Position);
        ax2 = axes(tl);
        ax2.XAxisLocation = 'top';
        ax2.YAxisLocation = 'right';
        ax2.YTickLabel = {};
        ax2.XMinorTick = 'on';
        ax2.Color = 'none';
        ax2.XLim = [OptArgs.Xmin, OptArgs.Xmax];
        ax2.YLim = ax.YLim;
        ax2.YScale = 'log';
        ax2.Box = 'off';
        xlabel(ax2,'Magnitide, $m$','Interpreter','latex');
    end
end
