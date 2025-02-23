%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 21 May 2022
%   ...
%   version 1.1.0, 26 August 2024
%

% extracting earthquakes from the catalogue
vCat0 = load_seismicity(Model.sEqCatName,'MagMin',Model.fMmin-0.5,'MagMax',Model.fMmax,'DateStart',Model.vDateStart,'Tstart',Model.fT0,'Tend',Model.fT1,...
                         'DepthMax',Model.fDepthMax,'SeismRegion',Model.vSeismRegion_targ);
vCat0(:,4) = round(vCat0(:,4)/Model.fDm)*Model.fDm;  % round the catalog to fDm however this eliminates the effect of 0.5*fDm when using fMmin
vCat = vCat0(vCat0(:,4)>=Model.fMmin,:); % this ensures that events with 0.5*fDm less than fMc are included
vTM = vCat(:,[1,4]);
