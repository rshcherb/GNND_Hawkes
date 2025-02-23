function mag_complete(vTM,dm,Model,varargin)
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 24 June 2022
%   ...
%   version: 1.0.0, 24 June 2022
%
    bExp          = false; vGOF_Exp = [];
    bPareto       = false; vGOF_Pareto = [];
    bTrnPareto    = false; vTrnPareto = [];
    bTapPareto    = false; vGOF_TapPareto = [];
    bTapTrnPareto = false; vTapTrnPareto = [];
    bTapTapPareto = false; vTapTapPareto = [];
    bTapParetoPareto = false; vTapParetoPareto = [];
    bParetoPareto = false; vParetoPareto = [];
    fAlpha = 0.05;
    bSaveFigure = false;
    for k = 1:length(varargin)
        if strcmp('Exp',varargin{k})
            bExp = true;
        end
        if strcmp('Pareto',varargin{k})
            bPareto = true;
        end
        if strcmp('TrnPareto',varargin{k})
            bTrnPareto = true;
        end
        if strcmp('TapPareto',varargin{k})
            bTapPareto = true;
        end
        if strcmp('TapTrnPareto',varargin{k})
            bTapTrnPareto = true;
        end
        if strcmp('TapTapPareto',varargin{k})
            bTapTapPareto = true;
        end
        if strcmp('TapParetoPareto',varargin{k})
            bTapParetoPareto = true;
        end
        if strcmp('ParetoPareto',varargin{k})
            bParetoPareto = true;
        end
        if strcmp('Alpha',varargin{k})
            fAlpha = varargin{k+1};
        end
        if strcmp('SaveFigure',varargin{k})
            %sFileName = varargin{k+1};
            bSaveFigure = true;
        end
    end
    
    vMag    = vTM(Model.nJs:Model.nJe,2);
    vMom    = mag2moment(vMag);
    Mom_min = mag2moment(Model.fMc);
    %vMagBin = linspace(Model.fMc,Model.m_max,Model.nBin)';
    vMagBin = (Model.fMc:dm:Model.m_max)';
    nBin = length(vMagBin);
    vMomBin = mag2moment(vMagBin);

    if bExp
        nPar = 2;
        vGOF_Exp = zeros(nBin,3);
        vGOF = zeros(nBin,nPar);
        for n = 1:nBin
            vMag_n = vMag(vMag >= vMagBin(n));
            [vPar, vErr] = gr_fit(vMag_n,vMagBin(n),Model.fMmax,Model.fDm,'Alpha',fAlpha);
            vPar(1) = log(10)*vPar(1); % convert from b-value to \beta
            gof = gof_ww(vMag_n,vMagBin(n),Model.fDm,'Exp',vPar(1));
            vGOF_Exp(n,:) = [vPar(1), vErr(1), gof];
            vGOF(n,:) = [gof, vPar(1)];
        end
        sModel = 'Exp';
        cYlabel = {'gof, \%','$\beta$'};
    end
    if bPareto
        nPar = 2;
        vGOF_Pareto = zeros(nBin,3);
        vGOF = zeros(nBin,nPar);
        for n = 1:nBin
            vMom_n = vMom(vMom >= vMomBin(n));
            vMag_n = vMag(vMag >= vMagBin(n));
            [vPar, vErr] = pareto_fit(vMom_n,vMomBin(n),'Alpha',fAlpha);
            vPar = log(10)*3/2*vPar;
            gof = gof_ww(vMag_n,vMagBin(n),Model.fDm,'Exp',vPar);
            vGOF_Pareto(n,:) = [vPar, vErr, gof];
            vGOF(n,:) = [gof, vPar];
        end
        sModel = 'Pareto';
        cYlabel = {'gof, \%','$\beta$'};
    end
    if bTrnPareto
        [vTrnPareto, vTrnParetoErr] = trn_pareto_fit(vMom,Mom_min,'Alpha',fAlpha);
    end
    if bTapPareto
        nPar = 3;
        vGOF_TapPareto = zeros(nBin,5);
        vGOF = zeros(nBin,nPar);
        for n = 1:nBin
            vMom_n = vMom(vMom >= vMomBin(n));
            %vMag_n = vMag(vMag >= vMagBin(n));
            [vPar, vErr] = tap_pareto_fit(vMom_n,vMomBin(n),'Alpha',fAlpha);
            gof = gof_ww(vMom_n,vMomBin(n),[],'TapPareto',vPar);
            vGOF_TapPareto(n,:) = [vPar, vErr, gof];
            vGOF(n,:) = [gof, vPar(1), moment2mag(vPar(2))];
        end
        sModel = 'Tap. Pareto';
        cYlabel = {'gof, \%','$\alpha$','$M_c$'};
    end
    if bTapTrnPareto
        gof = 0;
        vTapTrnPareto = [1.0 -0.61 1.34 2.86 0.91 0.09 gof]; % [alpha, m_cm, gamma, m_max, w1, w2, gof]
    end
    if bTapTapPareto
        [vTapTapPareto, vTapTapParetoErr] = taptap_pareto_fit(vMom,Mom_min,'Alpha',fAlpha);
    end
    if bTapParetoPareto
        [vTapParetoPareto, vTapParetoParetoErr] = tap_pareto_pareto_fit(vMom,Mom_min,'Alpha',fAlpha);
        vGOF_TapParetoPareto = zeros(Model.nBin,9);
        for n = 1:Model.nBin
            vMom_n = vMom(vMom >= vMomBin(n));
            %vMag_n = vMag(vMag >= vMagBin(n));
            [vPar, vErr] = tap_pareto_pareto_fit(vMom_n,vMomBin(n),'Alpha',fAlpha);
            gof = 0;
            vGOF_TapParetoPareto(n,:) = [vPar, vErr, gof];
        end
    end
    if bParetoPareto
        [vParetoPareto, vParetoParetoErr] = pareto_pareto_fit(vMom,Mom_min,'Alpha',fAlpha);
    end
    
    dm = 0.05;
    figure('Name','FMD completeness analysis');
    set(gcf,'color','w'); % background color for the figure

    for np = 1:nPar
        hsp{np} = subplot(nPar,1,np);
        plot(vMagBin,vGOF(:,np),'.','MarkerSize',8);
        ax = gca; ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
        grid on
        hold on
        ylabel(cYlabel{np},'Interpreter','latex');
    end
    plot(hsp{1},vMagBin,90*ones(nBin,1),'--');
    plot(hsp{1},vMagBin,95*ones(nBin,1),'-.');
    hold off
    sgtitle([sModel,': ',Model.sTitle],'FontSize',12);

%     if bExp
%         subplot(2,1,1)
%         plot(vMagBin,vGOF_Exp(:,1),'.','MarkerSize',8);
%         ax = gca; ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
%         grid on
%         ylabel('$\beta$','Interpreter','latex');
%         subplot(2,1,2)
%         plot(vMagBin,vGOF_Exp(:,3),'.','MarkerSize',8);
%         ax = gca; ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
%         grid on
%         hold on
%         plot(vMagBin,90*ones(nBin,1),'--');
%         plot(vMagBin,95*ones(nBin,1),'-.');
%         ylabel('gof, \%','Interpreter','latex');
%         hold off
%         sgtitle(['Exp: ',Model.sTitle],'FontSize',12);
%     end
%     if bPareto
%         subplot(2,1,1)
%         plot(vMagBin,vGOF_Pareto(:,1),'.');
%         ax = gca; ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
%         grid on
%         ylabel('$\beta$','Interpreter','latex');
%         subplot(2,1,2)
%         plot(vMagBin,vGOF_Pareto(:,3),'.');
%         ax = gca; ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
%         grid on
%         hold on
%         plot(vMagBin,90*ones(nBin,1),'--');
%         plot(vMagBin,95*ones(nBin,1),'-.');
%         ylabel('gof, \%','Interpreter','latex');
%         hold off
%         sgtitle(['Pareto: ',Model.sTitle],'FontSize',12);
%     end
%     if bTapPareto
%         subplot(3,1,1)
%         plot(vMagBin,vGOF_TapPareto(:,1),'.');
%         ax = gca; ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
%         grid on
%         ylabel('$\alpha$','Interpreter','latex');
%         subplot(3,1,2)
%         plot(vMagBin,moment2mag(vGOF_TapPareto(:,2)),'.');
%         ax = gca; ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
%         grid on
%         ylabel('$M_c$','Interpreter','latex');
%         subplot(3,1,3)
%         plot(vMagBin,vGOF_TapPareto(:,5),'.');
%         %ylim([-Inf,2]);
%         ax = gca; ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
%         grid on
%         hold on
%         plot(vMagBin,90*ones(nBin,1),'--');
%         plot(vMagBin,95*ones(nBin,1),'-.');
%         ylabel('gof, \%','Interpreter','latex');
%         hold off
%         sgtitle(['Tap. Pareto: ',Model.sTitle],'FontSize',12);
%     end
%     if bTapParetoPareto
%         subplot(4,1,1)
%         plot(vMagBin,vGOF_TapParetoPareto(:,1),'.');
%         ax = gca; ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
%         grid on
%         ylabel('$\alpha$','Interpreter','latex');
%         subplot(4,1,2)
%         plot(vMagBin,moment2mag(vGOF_TapParetoPareto(:,2)),'.');
%         ax = gca; ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
%         grid on
%         ylabel('$m_c$','Interpreter','latex');
%         subplot(4,1,3)
%         plot(vMagBin,vGOF_TapParetoPareto(:,3),'.');
%         ax = gca; ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
%         grid on
%         ylabel('$\beta$','Interpreter','latex');
%         subplot(4,1,4)
%         plot(vMagBin,vGOF_TapParetoPareto(:,4),'.');
%         ax = gca; ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
%         grid on
%         ylabel('$w$','Interpreter','latex');
%         sgtitle(['Tap. Pareto+Pareto: ',Model.sTitle],'FontSize',12);
%     end
    xlabel('Magnitude, $m_\mathrm{c}$','Interpreter','latex');

    if bSaveFigure
        sName = [Model.sFileName,'_mag_completeness'];
        saveas(gcf,[sName,'.png']);
        saveas(gcf,[sName,'.eps'],'epsc');
        %saveas(gcf,[sName,'.pdf']);
    end
end

