function prod = product_hawkes(vCat,fTs,fT,nJs,vPar)
%
%   Computes the productivity of the Hawkes process in the interval [fTs, fT] 
%   vCat  - events catalogue between fT0 and fTe 
%   fTs  - the start of the target time interval
%   fT   - the end of the time interval fT > fTs and less than or equal fTe
%   nJs  - the index of the first event in [fTs, fTe]
%   vPar - Hawkes parameters: [mu, A, alpha]
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 1 November 2024
%   ...
%   version: 1.0.0, 1 November 2024
%
    prod = sum(exp(-vPar(3).*(fTs - vCat(1:nJs-1,1))) - exp(-vPar(3).*(fT - vCat(1:nJs-1,1))));
    nJ = find(vCat(:,1) <= fT,1,'last'); % find the last index for which vCat(:,1) <= fT
    prod = vPar(1).*(fT - fTs) + vPar(2)./vPar(3).*(prod + sum( 1 - exp(-vPar(3).*(fT - vCat(nJs:nJ,1))) ) );
end

