%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 21 May 2022
%   ...
%   version 2.0.0, 31 October 2024
%
clear all
addpath(genpath('MATLAB/'));

sPointProc = 'Hawkes';   % 'ETAS';     % 
Model.fMc = -0.9;        % lower magnitude cutoff to perform analysis
common_parameters;
Model.sFileName = [sDir_rate,sFileName];
load_data;
if strcmp(sPointProc,'ETAS')
    Model = set_etas_model(vTM,Model,PPip,'OptimMethod',sSolver);
    Model.sTitle = [Model.sTitle0,': T_s = ', num2str(Model.fTs), ', T_e = ', num2str(Model.fTe),'; N = ', num2str(Model.nNum), '; m_c = ',num2str(Model.fMc),'; m_0 = ', num2str(Model.fM0)];
    Model.ETAS = model_etas_rate(vTM(1:Model.nJe,:),Model.ETAS,'PointProcess',sPointProc,'ErrorMethod',Model.sErrorMethod, ...
             'OptimMethod',sSolver,'Title',Model.sTitle,'SaveFigure',Model.sFileName);
elseif strcmp(sPointProc,'Hawkes')
    Model = set_hawkes_model(vTM,Model,PPip,'OptimMethod',sSolver);
    Model.sTitle = [Model.sTitle0,': T_s = ', num2str(Model.fTs), ', T_e = ', num2str(Model.fTe),'; N = ', num2str(Model.nNum), '; m_c = ',num2str(Model.fMc)];
    Model.Hawkes = model_hawkes_rate(vTM,Model.Hawkes,'PointProcess',sPointProc,'ErrorMethod',Model.sErrorMethod,...
                                  'OptimMethod',sSolver,'DateFormat',sDateFormat,'Title',Model.sTitle,'SaveFigure',Model.sFileName);
end
% frequncy-magnitude statistics for earthquakes used to fit the PointProc model
model_freq_mag(vTM(Model.nJs:Model.nJe,:),Model,Model.sTitle);
