function fGof = pp_gof(vCat,PP,gif_pp,OptArgs)
%
%   vCat   - the catalogue of earthquakes
%   PP     - input structure: nNumEq, fT0, fTs, fTe, fT1, fM0, vMs, ...
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 31 October 2024
%   ...
%   version: 1.0.0, 1 November 2024
%
    arguments
        vCat double
        PP struct                             % structure that defines the poit-process to fit
        gif_pp                                % function handle for the point process ground intensity function
        OptArgs.GoFmethod char = 'area'       % the method used to compute goodness-of-fit
        OptArgs.PointProcess char = 'Hawkes'  % 
    end

    vT        = vCat(PP.nJs:PP.nJe,1); % extract event times in the target window only
    vCumModel = gif_pp(vT,vCat,PP.vPar,1);
    indxNaN   = ~isnan(vCumModel); % it is possible to have NaN values when doing interpolation
    vCumModel = vCumModel(indxNaN);
    vCumData  = linspace(vCumModel(1),vCumModel(end),length(vCumModel))';
    if strcmp(OptArgs.GoFmethod,'area')
        fNormArea = 0.5*(vCumData(end) - vCumData(1))^2;   % area of the right triangle under the dioganal
        fGof  = trapz(vCumData,abs(vCumModel - vCumData))./fNormArea; % the normalized area of the difference between the cumulative curve in the transformed time and the dioganal 
    else
        fGoF  = 0;
    end
end

