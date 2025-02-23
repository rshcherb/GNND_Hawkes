%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 21 May 2022
%   ...
%   version 1.4.0, 4 November 2024
%
sDir_data         = 'data/';
Model.sEqCatName  = [sDir_data,'catalog.dat']; % catalog with blasts removed
Model.sSeismName  = 'Mine';
sDateFormat       = 'date';
Model.vSeismRegion_bg   = [2 lat1 lat2 lon1 lon2];
Model.vSeismRegion_targ = [2 lat1 lat2 lon1 lon2];

%   T0 ------- fTs ---------------- fTe ------- fT1
Model.vDateStart = [2022,1,1,0,0,0];        % 
Model.vDateTs    = [2022,3,1,0,0,0.0];      % 
Model.vDateEnd   = [2023,11,15,0,0,0];      % days since vDateStart

Model.fT0         = 0.0;                     % starting time
Model.fTs         = datenum(Model.vDateTs) - datenum(Model.vDateStart);  % start of the target time interval
Model.fTe         = datenum(Model.vDateEnd) - datenum(Model.vDateStart); % end of the target time interval
Model.fDT         = 1.0;                    % forecasting time interval in days  
Model.sRateNorm   = 'Harte';      % 'Ogata'; %  
Model.sErrorMethod = 'Hessian';   % 'Simul'; % 'BootStrap'; % 'LogLik';    % 'OVP'; % 
sSolver           = 'fmincon';    % 'fminsearch'; % 'gs'; % 
if ~exist('sPointProc','var')
    sPointProc    = 'ETAS';       % point process to use
end
PPip.Aeq = []; PPip.beq = [];     % no constraints on parameters
if strcmp(sPointProc,'ETAS')
    Model.sNamePart = '_etas';
    % set the initial parameters for the ETAS fitting
    Model.vParName   = {'\mu'  'K'   'c'    'p'   '\alpha'};
    PPip.vInitParams = [0.1     0.1   0.5    1.1   0.6];
    PPip.vLowB    =    [0.0     0.01  0.0    0.0   0.0];
    PPip.vUppB    =    [100.0   100.0  5.0    5.0   3.0];
%     PPip.Aeq = zeros(5); PPip.Aeq(1,1) = 1; %
%     PPip.beq = zeros(5,1); PPip.beq(1) = 0.0; % mu
    sDir_rate     = 'eq_rate_etas/';
elseif strcmp(sPointProc,'Hawkes')
    Model.sNamePart = '_hawkes';
    % set the initial parameters for the Hawkes fitting
    Model.vParName   = {'\mu'  'A'   '\alpha'};
    PPip.vInitParams = [0.5     0.1   0.2];
    PPip.vLowB    =    [0.0     0.0   0.0];
    PPip.vUppB    =    [1e2     1e2   1e2];
%     PPip.Aeq = zeros(3); PPip.Aeq(1,1) = 1; %
%     PPip.beq = zeros(3,1); PPip.beq(1) = 0.0; % mu
    sDir_rate     = 'eq_rate_hawkes/';
end
Model.nPar        = 1 + length(PPip.vInitParams); % the total number of parameters beta + ETAS parameters
Model.fT1         = Model.fTe + Model.fDT;
Model.fDm         = 0.01;
Model.fAlpha      = 0.05;
Model.fMmin       = Model.fMc - 0.5*Model.fDm;  % lower cutoff for magnitudes to extract or generate
Model.fMmax       = 2.3;                        % upper cutoff for magnitudes to use
Model.fM0         = Model.fMc;                  % reference magnitude for the ETAS model
Model.fDepthMin   = -15.0;
Model.fDepthMax   = 15.0;

sDateRange        = [num2str(Model.vDateStart(1)),'/',num2str(Model.vDateStart(2),2),'/',num2str(Model.vDateStart(3),2),'-',...
                     num2str(Model.vDateEnd(1)),'/',num2str(Model.vDateEnd(2),2),'/',num2str(Model.vDateStart(3),2)];
Model.sTitle0     = [Model.sSeismName,': ',sDateRange];
sFileName         = [Model.sSeismName,'_',replace(sDateRange,"/","_"),'_',num2str(Model.fMc),'m_',num2str(Model.fTs),'d_',num2str(Model.fTe),'d'];  % output file name


