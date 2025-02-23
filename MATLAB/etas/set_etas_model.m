function Model = set_etas_model(vCat,Model,mIP,OptArgs)
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 6 May 2021
%   ...
%   version 2.0.0, 3 November 2024
%
    arguments
        vCat double
        Model struct
        mIP struct
        OptArgs.PointProcess char = 'ETAS'  % 
        OptArgs.OptimMethod char = 'fmincon'       % 'gs'; 'fminsearch';  
    end
    
    Model.ETAS = set_model2pp(vCat,Model,mIP,'OptimMethod',OptArgs.OptimMethod);
    
    Model.nJs     = Model.ETAS.nJs;
    Model.nJe     = Model.ETAS.nJe;
    %Model.nNum    = length(vTM(Model.nJs:Model.nJe,1)); % number of earthquakes between fTs and fTe
    Model.nNum    = Model.ETAS.nNum; % number of earthquakes between fTs and fTe

    if ~isfield(Model,'sPointProc')
        Model.sPointProc = 'ETAS';  % the default point process
    end
    if strcmp(Model.sPointProc,'ETASFrackConv') || strcmp(Model.sPointProc,'ETAS_MultConv')
        if strcmp(Model.sConvKernel,'Pareto')
            Model.ETAS.nConvKernel = 1;
        elseif strcmp(Model.sConvKernel,'ExpPdf')
            Model.ETAS.nConvKernel = 2;
        elseif strcmp(Model.sConvKernel,'PowerLaw')
            Model.ETAS.nConvKernel = 3;
        end
        Model.ETAS.sConvKernel = Model.sConvKernel;
    end
end
