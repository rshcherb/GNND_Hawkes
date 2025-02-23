function vEQseq = load_seismicity(sEqCatName,varargin)
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 8 January 2021
%   ...
%   version 1.0.0, 8 January 2021
%
    fMmin         = 0.0;   % can be negative 
    fMmax         = 9.9;
    fDepthMin     = -10.0; % depath is negative upwards
    fDepthMax     = 700.0; % depth is positive downwards
    vDateStart    = [1900 1 1 0.0];
    vDateEnd      = [2100 1 1 0.0];
    fT0           = 0.0;
    fT1           = 100.0;
    vSeismRegion  = [2 -90 90 -360 360];
    bCatFileOut2  = false;
    bCatFileOut8  = false;
    bCatFileOut10 = false;
    for k = 1:length(varargin)
        if strcmp('MagMin',varargin{k})
            fMmin = varargin{k+1};
        end
        if strcmp('MagMax',varargin{k})
            fMmax = varargin{k+1};
        end
        if strcmp('DateStart',varargin{k})
            vDateStart = varargin{k+1};
        end
        if strcmp('DateEnd',varargin{k})
            vDateEnd = varargin{k+1};
        end
        if strcmp('Tstart',varargin{k})
            fT0 = varargin{k+1};
        end
        if strcmp('Tend',varargin{k})
            fT1 = varargin{k+1};
        end
        if strcmp('SeismRegion',varargin{k})
            vSeismRegion = varargin{k+1};
        end
        if strcmp('DepthMin',varargin{k})
            fDepthMin = varargin{k+1};
        end
        if strcmp('DepthMax',varargin{k})
            fDepthMax = varargin{k+1};
        end
        if strcmp('CatFileOut2',varargin{k})
            bCatFileOut2 = true;
            sCatFileOut = varargin{k+1};
        end
        if strcmp('CatFileOut8',varargin{k})
            bCatFileOut8 = true;
            sCatFileOut = varargin{k+1};
        end
        if strcmp('CatFileOut10',varargin{k})
            bCatFileOut10 = true;
            sCatFileOut = varargin{k+1};
        end
    end
    %disp(fMmin)
    rawcat = load(sEqCatName,'-ascii');   % loading a user specified catalog file with 10 columns
    indx = (rawcat(:,9) >= fMmin) & (rawcat(:,9) <= fMmax) & (rawcat(:,10) <= fDepthMax);
    mEQcat(:,1) = datenum(rawcat(indx,1:6));
    mEQcat(:,2:5) = rawcat(indx,7:10);
    
    pd          = datenum(vDateStart);    % the pivot time given by vDateStart
    fTimeMin    = pd + fT0;               % the reference date corresponding to T0
    fTimeMax    = pd + fT1;               % the end of the target time interval
    indx        = (mEQcat(:,1) >= fTimeMin)  & (mEQcat(:,1) <= fTimeMax); % indx for earthquakes between fTimeMin and fTimeMax
    mEQcat      = mEQcat(indx,:);
    vEQseq      = seism_in_georegion(mEQcat,fDepthMin,fDepthMax,vSeismRegion);
    vEQseq(:,1) = vEQseq(:,1) - pd;       % subtract the pivot time to make all earthquakes times relative to pd
    
    if bCatFileOut2
        cindx = [1, 0, 0, 1, 0] == 1; % select time and magnitude columns
        vTM = vEQseq(:,cindx);
        save([sCatFileOut,'_eq_cat_2.dat'],'vTM','-ascii');
    end
    if bCatFileOut8
        fid = fopen([sCatFileOut,'_eq_cat_8.dat'],'w');
        for n = 1:length(vEQseq(:,1))
            dv = datevec(pd + vEQseq(n,1));
            year  = dv(1);
            month = dv(2);
            day   = dv(3);
            hour  = dv(4) + dv(5)/60 + dv(6)/3600;
            lat = vEQseq(n,2);
            lon = vEQseq(n,3);
            mag = vEQseq(n,4);
            dep = vEQseq(n,5);
            fprintf(fid,'%4d %2d %2d %9.6f %8.4f %9.4f %4.2f %6.2f\n',...
                year,month,day,hour,lat,lon,mag,dep);
        end
        fclose(fid);
    end
    if bCatFileOut10
        fid = fopen([sCatFileOut,'_eq_cat_10.dat'],'w');
        for n = 1:length(vEQseq(:,1))
            dv = datevec(pd + vEQseq(n,1));
            year  = dv(1);
            month = dv(2);
            day   = dv(3);
            hour  = dv(4);
            min   = dv(5);
            sec   = dv(6);
            lat = vEQseq(n,2);
            lon = vEQseq(n,3);
            mag = vEQseq(n,4);
            dep = vEQseq(n,5);
            fprintf(fid,'%4d %2d %2d %2d %2d %5.2f %8.4f %9.4f %4.2f %6.2f\n',...
                year,month,day,hour,min,sec,lat,lon,mag,dep);
        end
        fclose(fid);
    end
end

