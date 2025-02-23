%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 21 May 2022
%   ...
%   version 1.0.1, 28 August 2024
%
clear all
addpath(genpath('MATLAB/'));

Model.fMc = -1.0;         % lower magnitude cutoff to perform analysis
common_parameters;
sDir_out = 'fm/';
Model.sFileName = [sDir_out,sFileName];
load_data;
Model = set_etas_model(vTM,Model,PPip,'OptimMethod',sSolver);
Model.sTitle = [Model.sTitle0,': T_s = ', num2str(Model.fTs), ', T_e = ', num2str(Model.fTe),'; N = ', num2str(Model.nNum), '; m_c = ',num2str(Model.fMc)];

% frequncy-magnitude statistics for earthquakes used to fit the PointProc model
model_freq_mag(vTM(Model.nJs:Model.nJe,:),Model,Model.sTitle);

%model_FMD(vTM,Model,'Pareto','TapPareto','TapTapPareto','TapParetoPareto','SaveFigure');
model_FMD(vTM,Model,'Pareto','TapTapPareto','SaveFigure');

%model_FMD(vTM,Model,'Pareto','TapPareto','TapTapPareto','TapParetoPareto','ParetoPareto','SaveFigure');
%model_FMD(vTM,Model,'Pareto','TapTapPareto','TapParetoPareto','SaveFigure');
