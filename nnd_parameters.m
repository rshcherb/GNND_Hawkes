%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 17 September 2024
%   ...
%   version 1.2.0, 27 October 2024
%
NNDpar.MixtModel = 'GMM';
NNDpar.ThreshMethod = 'intersect';
NNDpar.sEtaDistrbFit = 'Normal'; % the model to fit to each mode of \eta distribution: 'Weibull', 'Normal'
NNDpar.b = 1.0;
NNDpar.df = 1.6;
NNDpar.q = 0.5;
NNDpar.r_min = 0.0;
NNDpar.r_max = 20.0;     % maximum radial distance in km. This is used in nnd() function
NNDpar.dt_max = 100;     % maximum time in days betweet parent-daughter pairs. This is used for declustering in nnd_decluster_eta()
NNDpar.dr_max = 20.0;     % maximum distance in km betweet parent-daughter pairs
NNDpar.nGMMmax = 2;
NNDpar.nNode_min = 10;    % the minimum number of nodes above which to plot
NNDpar.nClst_largest = 9; % the number of the largest clusters to plot
if strcmp (sFMD,'GR')
    NNDpar.EtaLim = [10^(-5.0), 10^(2)];
    NNDpar.EtaYLim = [0, 0.7];
    NNDpar.TLim = [10^(-4), 10^(4)];
    NNDpar.RLim = [10^(-4), 10^(3)];
elseif strcmp(sFMD,'Exp')
    NNDpar.EtaLim = [10^(-5.0), 10^(2)];
    NNDpar.EtaYLim = [0, 0.7];
    NNDpar.TLim = [10^(-4.5), 10^(4)];
    NNDpar.RLim = [10^(-5.0), 10^(3.5)];
    NNDpar.beta = 1.0*log(10);
    NNDpar.m_min = Model.fMc;
elseif strcmp(sFMD,'TapTapPareto')
    NNDpar.EtaLim = [10^(-5.0), 10^(2)];
    NNDpar.EtaYLim = [0, 0.7];
    NNDpar.TLim = [10^(-4.5), 10^(3)];
    NNDpar.RLim = [10^(-5.0), 10^(3.5)];
elseif strcmp(sFMD,'TapParetoPareto')
    NNDpar.EtaLim = [10^(-5.0), 10^(2)];
    NNDpar.EtaYLim = [0, 0.7];
    NNDpar.TLim = [10^(-4.5), 10^(3)];
    NNDpar.RLim = [10^(-5.0), 10^(3.5)];
end
% duplicate in Model
Model.sEtaDistrbFit = NNDpar.sEtaDistrbFit; % the model to fit to each mode of \eta distribution: 'Weibull', 'Normal'
%Model.EtaLim = NNDpar.EtaLim;
%Model.EtaYLim = NNDpar.EtaYLim;
%Model.TLim = NNDpar.TLim;
%Model.RLim = NNDpar.RLim;
Model.nGMMmax = NNDpar.nGMMmax;
Model.nNode_min = NNDpar.nNode_min;    % the minimum number of nodes above which to plot
Model.nClst_largest = NNDpar.nClst_largest; % the number of the largest clusters to plot

%Model.NNDpar = NNDpar;

