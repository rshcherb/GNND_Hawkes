function model_FMD(vTM,Model,varargin)
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 18 June 2022
%   ...
%   version: 2.0.1, 16 November 2024
%
    bPareto       = false; Pareto = [];
    bTrnPareto    = false; TrnPareto = [];
    bTapPareto    = false; TapPareto = [];
    bTapTrnPareto = false; TapTrnPareto = [];
    bTapTapPareto = false; TapTapPareto = [];
    bTapParetoPareto = false; TapParetoPareto = [];
    bParetoPareto = false; ParetoPareto = [];
    fAlpha = 0.05;
    bSaveFigure = false;
    bDisplay = false;
    for k = 1:length(varargin)
        if strcmp('Pareto',varargin{k}),          bPareto = true; end
        if strcmp('TrnPareto',varargin{k}),       bTrnPareto = true; end
        if strcmp('TapPareto',varargin{k}),       bTapPareto = true; end
        if strcmp('TapTrnPareto',varargin{k}),    bTapTrnPareto = true; end
        if strcmp('TapTapPareto',varargin{k}),    bTapTapPareto = true; end
        if strcmp('TapParetoPareto',varargin{k}), bTapParetoPareto = true; end
        if strcmp('ParetoPareto',varargin{k}),    bParetoPareto = true; end
        if strcmp('Alpha',varargin{k}), fAlpha = varargin{k+1}; end
        if strcmp('SaveFigure',varargin{k})
            %sFileName = varargin{k+1};
            bSaveFigure = true;
        end
        if strcmp('Display',varargin{k}), bDisplay = varargin{k+1}; end
    end
    vPos = [400 300 800 600];

    vMag = vTM(Model.nJs:Model.nJe,2);
    vMom = mag2moment(vMag);
    fMom_min = mag2moment(Model.fMc);

    if bPareto
        [vPar, vPci, vPerr] = pareto_fit(vMom,fMom_min,'Alpha',fAlpha);
        if bDisplay
            fprintf('%s\n',sprintf('Pareto: beta = %.3f \xB1 %.3f',vPar(1),vPerr(1)));
        end
        gof = gof_ww(vMom,fMom_min,'Pareto',vPar);
        %gof = gof_ww(vMag,Model.fMc,'Exp',1.5*vPar(1)*log(10.0),'Units','magnitude'); % uisng the magnitude distribution and exponential function
        Pareto.par = vPar;
        Pareto.pci = vPci;
        Pareto.err = vPerr;
        Pareto.gof = gof;
        nPar = 1;
        vPdf = pareto_pdf(vMom,vPar(1),fMom_min);
        fLLe = sum(log(vPdf)); % using moment
        Pareto.AICmom = -2.0*fLLe + 2*nPar;
        fLLe = sum(log(pdf_moment2mag(vMag,vPdf))); % using magnitude
        Pareto.AICmag = -2.0*fLLe + 2*nPar;
    end
    if bTrnPareto
        [vPar, vPci, vPerr] = trn_pareto_fit(vMom,fMom_min,'Alpha',fAlpha);
        magerr = moment2mag(vPar(2) + vPerr(2)) - moment2mag(vPar(2) - vPerr(2));
        if bDisplay
            fprintf('%s\n',sprintf('Truncated Pareto: gamma = %.3f \xB1 %.3f',vPar(1),vPerr(1)))
            fprintf('%s\n',sprintf('                  M_max = %.3e \xB1 %.3e',vPar(2),vPerr(2)))
            fprintf('%s\n',sprintf('                  m_max = %.3f \xB1 %.3f',moment2mag(vPar(2)),magerr))
        end
        TrnPareto.par = vPar;
        TrnPareto.pci = vPci;
        TrnPareto.err = vPerr;
        TrnPareto.magerr = magerr;
        TrnPareto.gof = 0;
        nPar = 2;
        vPdf = trn_pareto_pdf(vMom,vPar(1),fMom_min,vPar(2));
        fLLe = sum(log(vPdf)); % using moment
        TrnPareto.AICmom = -2.0*fLLe + 2*nPar;
        fLLe = sum(log(pdf_moment2mag(vMag,vPdf))); % using magnitude
        TrnPareto.AICmag = -2.0*fLLe + 2*nPar;
    end
    if bTapPareto
        [vPar, vPci, vPerr] = tap_pareto_fit(vMom,fMom_min,'Alpha',fAlpha);
        mcm_err = moment2mag(vPar(2) + vPerr(2)) - moment2mag(vPar(2) - vPerr(2));
        if bDisplay
            fprintf('%s\n',sprintf('Tapered Pareto: alpha = %.3f \xB1 %.3f',vPar(1),vPerr(1)))
            fprintf('%s\n',sprintf('                M_cm = %.3e \xB1 %.3e',vPar(2),vPerr(2)))
            fprintf('%s\n',sprintf('                m_cm = %.3f \xB1 %.3f',moment2mag(vPar(2)),mcm_err))
        end
        TapPareto.par = vPar;
        TapPareto.pci = vPci;
        TapPareto.err = vPerr;
        TapPareto.magerr = mcm_err;
        TapPareto.gof = gof_ww(vMom,fMom_min,'TapPareto',vPar);
        nPar = 2;
        vPdf = tap_pareto_pdf(vMom,vPar(1),fMom_min,vPar(2));
        fLLe = sum(log(vPdf)); % using moment
        TapPareto.AICmom = -2.0*fLLe + 2*nPar;
        fLLe = sum(log(pdf_moment2mag(vMag,vPdf))); % using magnitude
        TapPareto.AICmag = -2.0*fLLe + 2*nPar;
    end
    if bTapTrnPareto
        vPar = [1.0 -0.61 1.34 2.86 0.91 0.09]; % [alpha, m_cm, gamma, m_max, w1, w2]
        vPerr = zeros(size(vPar));
        magerr = moment2mag(vPerr(2));
        TapTrnPareto.par = vPar;
        TapTrnPareto.err = vPerr;
        TapTrnPareto.magerr = magerr;
        TapTrnPareto.gof = 0;
        TapTrnPareto.AIC = 0;
    end
    if bTapTapPareto
        [vPar, vPci, vPerr] = taptap_pareto_fit(vMom,fMom_min,'Alpha',fAlpha);
        mcm_err(1) = moment2mag(vPar(2) + vPerr(2)) - moment2mag(vPar(2) - vPerr(2));
        mcm_err(2) = moment2mag(vPar(4) + vPerr(4)) - moment2mag(vPar(4) - vPerr(4));
        if bDisplay
            fprintf('%s\n',sprintf('Tap.+Tap. Pareto: alpha_1 = %.3f \xB1 %.3f',vPar(1),vPerr(1)))
            fprintf('%s\n',sprintf('                  M_cm_1 = %.3e \xB1 %.3e',vPar(2),vPerr(2)))
            fprintf('%s\n',sprintf('                  m_cm_1 = %.3f \xB1 %.3f',moment2mag(vPar(2)),mcm_err(1)))
            fprintf('%s\n',sprintf('                  alpha_2 = %.3f \xB1 %.3f',vPar(3),vPerr(3)))
            fprintf('%s\n',sprintf('                  M_cm_2 = %.3e \xB1 %.3e',vPar(4),vPerr(4)))
            fprintf('%s\n',sprintf('                  m_cm_2 = %.3f \xB1 %.3f',moment2mag(vPar(4)),mcm_err(2)))
            fprintf('%s\n',sprintf('                  w = %.3f \xB1 %.3f',vPar(5),vPerr(5)))
        end
        TapTapPareto.par = vPar;
        TapTapPareto.err = vPerr;
        TapTapPareto.pci = vPci;
        TapTapPareto.magerr = mcm_err;
        TapTapPareto.gof = gof_ww(vMom,fMom_min,'TapTapPareto',vPar);
        %fLLe = taptap_pareto_llfun(vMom,Mom_min,vPar(1),vPar(2),vPar(3),vPar(4),vPar(5));
        vPdf = taptap_pareto_pdf(vMom,fMom_min,vPar(1),vPar(2),vPar(3),vPar(4),vPar(5));
        nPar = 5;
        fLLe = sum(log(vPdf)); % using moment
        TapTapPareto.AICmom = -2.0*fLLe + 2*nPar;
        fLLe = sum(log(pdf_moment2mag(vMag,vPdf))); % using magnitude
        TapTapPareto.AICmag = -2.0*fLLe + 2*nPar;
    end
    if bTapParetoPareto
        [vPar, vPci, vPerr] = tap_pareto_pareto_fit(vMom,fMom_min,'Alpha',fAlpha);
        mcm_err = moment2mag(vPar(2) + vPerr(2)) - moment2mag(vPar(2) - vPerr(2));
        if bDisplay
            fprintf('%s\n',sprintf('Tap. Pareto+Pareto: alpha = %.3f \xB1 %.3f',vPar(1),vPerr(1)))
            fprintf('%s\n',sprintf('                      M_cm = %.3e \xB1 %.3e',vPar(2),vPerr(2)))
            fprintf('%s\n',sprintf('                      m_cm = %.3f \xB1 %.3f',moment2mag(vPar(2)),mcm_err))
            fprintf('%s\n',sprintf('                      beta = %.3f \xB1 %.3f',vPar(3),vPerr(3)))
            fprintf('%s\n',sprintf('                      w = %.3f \xB1 %.3f',vPar(4),vPerr(4)))
        end
        TapParetoPareto.par = vPar;
        TapParetoPareto.pci = vPci;
        TapParetoPareto.err = vPerr;
        TapParetoPareto.magerr = mcm_err;
        TapParetoPareto.gof = gof_ww(vMom,fMom_min,'TapParetoPareto',vPar);
        vPdf = tap_pareto_pareto_pdf(vMom,fMom_min,vPar(1),vPar(2),vPar(3),vPar(4));
        nPar = 4;
        fLLe = sum(log(vPdf)); % using moment
        TapParetoPareto.AICmom = -2.0*fLLe + 2*nPar;
        fLLe = sum(log(pdf_moment2mag(vMag,vPdf))); % using magnitude
        TapParetoPareto.AICmag = -2.0*fLLe + 2*nPar;
    end
    if bParetoPareto
        [vPar, vPci, vPerr] = pareto_pareto_fit(vMom,fMom_min,'Alpha',fAlpha);
        mcm_err = moment2mag(vPerr(2));
        if bDisplay
            fprintf('%s\n',sprintf('Pareto+Pareto: gamma = %.3f \xB1 %.3f',vPar(1),vPerr(1)))
            fprintf('%s\n',sprintf('               M_cm = %.3e \xB1 %.3e',vPar(2),vPerr(2)))
            fprintf('%s\n',sprintf('               m_cm = %.3f \xB1 %.3f',moment2mag(vPar(2)),mcm_err))
            fprintf('%s\n',sprintf('               beta = %.3f \xB1 %.3f',vPar(3),vPerr(3)))
        end
        ParetoPareto.par = vPar;
        ParetoPareto.pci = vPci;
        ParetoPareto.err = vPerr;
        ParetoPareto.magerr = mcm_err;
        ParetoPareto.gof = 0;
        ParetoPareto.AIC = 0;
    end
    
    figure('Name','Fits of the frequency-moment models','Position',vPos);
    % plot frequency-moment distributions
    [vX, vHist_emp, vPdf_emp, vCdf_emp] = plot_FMD(vMag,Model.fMc,Model.fDm,'Alpha',Model.fAlpha,'Units','moment',...
                                              'Pareto',Pareto,'TapPareto',TapPareto,'TrnPareto',TrnPareto,'TapTrnPareto',TapTrnPareto,...
                                              'TapTapPareto',TapTapPareto,'TapParetoPareto',TapParetoPareto,'ParetoPareto',ParetoPareto);
                                              %'Xmin',10^13,'Xmax',10^17);
    %title(Model.sTitle);
    if bSaveFigure
        sName = [Model.sFileName,'_moment_df'];
        save_cf(gcf,sName,'fig','png','pdf');
    end

    figure('Name','Fits of the frequency-magniude models','Position',vPos);
    % plot frequency-magnitude distributions
    [vX, vHist_emp, vPdf_emp, vCdf_emp] = plot_FMD(vMag,Model.fMc,Model.fDm,'Alpha',Model.fAlpha,'Units','magnitude',...
                                              'Pareto',Pareto,'TapPareto',TapPareto,'TrnPareto',TrnPareto,'TapTrnPareto',TapTrnPareto,...
                                              'TapTapPareto',TapTapPareto,'TapParetoPareto',TapParetoPareto,'ParetoPareto',ParetoPareto);
    if bSaveFigure
        sName = [Model.sFileName,'_magnitude_df'];
        save_cf(gcf,sName,'fig','png','pdf');
    end
end

