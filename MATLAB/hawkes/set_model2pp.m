function PP = set_model2pp(vCat,Model,mIP,OptArgs)
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 31 October 2024
%   ...
%   version 1.0.0, 31 October 2024
%
    arguments
        vCat double
        Model struct
        mIP struct
        OptArgs.OptimMethod char = 'fmincon'       % 'gs'; 'fminsearch';  
    end
    
%   T0 ------- Ts ---------------- Te ----- T1
    PP.vDateStart  = Model.vDateStart;
    PP.fT0         = Model.fT0;     % 
    PP.fTs         = Model.fTs;     % 
    PP.fTe         = Model.fTe;     %
    PP.fT1         = Model.fT1;     % 
    PP.fM0         = Model.fM0;     % the reference magnitude
    PP.fMc         = Model.fMc;     % the lower magnitude cutoff
    PP.fMmin       = Model.fMmin;   %
    PP.fMmax       = Model.fMmax;   %
    PP.mIP         = mIP;                 % the initial parameters
    PP.OptimMethod = OptArgs.OptimMethod; % optimization method
    PP.fAlpha      = Model.fAlpha;        % confidence level
    PP.nJs         = find(vCat(:,1) >= Model.fTs,1,'first'); % the first event whose time is greater than or equal to fTs
    PP.nJe         = find(vCat(:,1) <= Model.fTe,1,'last');  % the last event whose time is less than or equal fTe
    PP.nNtot       = length(vCat(:,1));                      % number of events in the catalogue vTM
    PP.nNum        = PP.nJe - PP.nJs + 1;                    % the number of earthquakes between fTs and fTe
    Model.fMagMean = mean(vCat(PP.nJs:PP.nJe,2));            % sample mean magnitude only for magnitudes between nJs and nJe
    PP.ParName     = Model.vParName;
end
