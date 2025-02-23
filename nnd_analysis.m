%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 25 May 2022
%   ...
%   version 1.3.0, 3 November 2024
%
clear all
addpath(genpath('MATLAB/'));

Model.fMc = -1.0;           % lower magnitude cutoff to perform analysis
sFMD      = 'TapTapPareto'; % 'TapParetoPareto'; % 'Exp';          % 'GR';           % frequency-magnitude distribution to rescale \eta, T, and R
common_parameters;
nnd_parameters;
sDir_out  = 'nnd/';
Model.sFileName = [sDir_out,sFileName,'_',sFMD];
load_data;
Model = set_etas_model(vCat,Model,PPip,'OptimMethod',sSolver);
Model.sTitle = [Model.sTitle0,': T_s = ', num2str(Model.fTs), ', T_e = ', num2str(Model.fTe),'; N = ', num2str(Model.nNum), '; m_c = ',num2str(Model.fMc)];
%Model.sTitle = [];

if strcmp(sFMD,'TapTapPareto') || strcmp(sFMD,'TapParetoPareto')
    vMag = vCat(Model.nJs:Model.nJe,4);
    vMom = mag2moment(vMag);
    fMom_min = mag2moment(Model.fMc);
    if strcmp(sFMD,'TapTapPareto')
        [vPar, vParErr] = taptap_pareto_fit(vMom,fMom_min,'Alpha',Model.fAlpha); % vPar = [M_min, alpha1, M_cm1, alpha2, M_cm2, w]
        NNDpar.M_min = fMom_min;
        NNDpar.alpha1 = vPar(1);
        NNDpar.M_cm1  = vPar(2);
        NNDpar.alpha2 = vPar(3);
        NNDpar.M_cm2  = vPar(4);
        NNDpar.w      = vPar(5);
    elseif strcmp(sFMD,'TapParetoPareto')
        [vPar, vParErr] = tap_pareto_pareto_fit(vMom,fMom_min,'Alpha',Model.fAlpha); % vPar = [M_min, alpha, M_cm, beta, w]
        NNDpar.M_min  = fMom_min;
        NNDpar.alpha  = vPar(1);
        NNDpar.M_cm   = vPar(2);
        NNDpar.beta   = vPar(3);
        NNDpar.w      = vPar(4);
    end
end
Model.NNDpar = NNDpar;

% nnd declustering
model_nnd(vCat,NNDpar,Model,'FMD',sFMD);

% frequency-magnitude statistics for earthquakes used 
model_freq_mag(vTM(Model.nJs:Model.nJe,:),Model,Model.sTitle);

model_FMD(vTM,Model,sFMD,'SaveFigure');
