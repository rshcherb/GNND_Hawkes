function vRate = gif_hawkes(t,vCat,fTs,nJs,vPar,bRateType)
%
%   Computes the conditional Hawkes rate at given times vT 
%   t   - times at which to compute the rate: fTs <= vT <= fTe
%   vCat  - events catalogue between fT0 and fT1. Assumes vTM(:,2) - fM0
%   fTs  - the start of the target time interval
%   nJs  - the index of the first event in [fTs, fTe]
%   vPar - Hawkes parameters: [mu, A, alpha]
%   bRateType - 0 - 'rate'; 1 - 'cumulative'
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 31 October 2024
%   ...
%   version: 1.0.0, 31 October 2024
%
    nt = length(t);
    vRate = zeros(nt,1);
    if bRateType == 0    % 0 - 'rate'
        for j = 1:nt
            nJ = find(vCat(:,1) < t(j),1,'last'); % find the last index for which vCat(:,1) < t(j)
            vRate(j) = vPar(1) + vPar(2).*sum( exp(-vPar(3).*(t(j) - vCat(1:nJ,1))) );
        end
    else                 % 1 - 'cumulative'
        for j = 1:nt
            vRate(j) = product_hawkes(vCat,fTs,t(j),nJs,vPar);
        end
    end
end

