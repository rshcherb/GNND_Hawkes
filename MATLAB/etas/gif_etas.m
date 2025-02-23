function vRate = gif_etas(vt,vCat,fTs,nJs,vPar,bRateType)
%
%   Computes the conditional ETAS rate at given times vT 
%   t   - times at which to compute the rate: fTs <= vT <= fTe
%   vCat  - events catalogue between fT0 and fT1. Assumes vTM(:,2) - fM0
%   fTs  - the start of the target time interval
%   nJs  - the index of the first event in [fTs, fTe]
%   vPar - Hawkes parameters: [mu, A, alpha]
%   bRateType - 0 - 'rate'; 1 - 'cumulative'
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 1 September, 2018
%   ...
%   version: 2.1.0, 2 November 2024
%
    nt = length(vt);
    vRate = zeros(nt,1);
    if bRateType == 0    % 0 - 'rate'
        for j = 1:nt
            nJ = find(vCat(:,1) < vt(j),1,'last'); % find the last index for which vCat(:,1) < t(j)
            vRate(j) = vPar(1) + vPar(2).*sum( exp(vPar(5).*vCat(1:nJ,2) )./( ((vt(j) - vCat(1:nJ,1))./vPar(3) + 1).^vPar(4) ) );
        end
    else                 % 1 - 'cumulative'
        for j = 1:nt
            vRate(j) = vPar(1).*(vt(j) - fTs) + product_etas(vCat,fTs,vt(j),nJs,vPar(2:5));
        end
    end
end

