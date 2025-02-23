function mEQseq = seism_in_georegion(mEQcat,fDepthMin,fDepthMax,vReg,varargin)
%
%   Finds all the earthquakes inside a given region
%   mEQcat - earthquake catalog
%   vReg   - the region
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 17 November 2019
%   ...
%   version 1.2.0, 21 November 2019
%
    sRegionType = 'Polygon';
    for k = 1:length(varargin)
        if strcmp('Region',varargin{k})
            sRegionType = varargin{k+1};
        end
    end

    indx = (mEQcat(:,5) >= fDepthMin) & (mEQcat(:,5) <= fDepthMax);
    mEQcat = mEQcat(indx,:);
    
    if strcmp(sRegionType,'Polygon')
        % compute the polygon corresponding to specific shapes
        [vLat, vLon] = georegion(vReg);
        % find points inside the polygon
        [indx_reg,on] = inpolygon(mEQcat(:,3),mEQcat(:,2),vLon,vLat);
    elseif strcmp(sRegionType,'Polyhedron')
        % one can use boundary() or alphaShape() to creat 3D region
        % then use https://www.mathworks.com/matlabcentral/fileexchange/37856-inpolyhedron-are-points-inside-a-triangulated-volume
        % or https://www.mathworks.com/matlabcentral/fileexchange/10226-inhull
        % or https://www.mathworks.com/matlabcentral/answers/101396-is-there-a-function-in-matlab-for-detecting-points-inside-a-polyhedron
        % using tsearchn()
    end
    % update the catalog
    mEQseq = mEQcat(indx_reg,:);
end

