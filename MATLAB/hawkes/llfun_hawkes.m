function fLle = llfun_hawkes(vPar,vCat,fTs,fTe,nJs,nJe)
%
%   The log-likelihood function for the Hawkes conditional rate
%
%   vPar     - the parameters of the Hawkes model as variables: \mu = vPar(1); A = vPar(2); \alpha = vPar(3);
%   vCat     - earthquake times and magnitudes with the subtracted reference magnitude m_0 (vTM(:,2) - m_0)
%   fTs      - the start time for the target window
%   fTe      - the end time for the target window
%   nJs      - the first index of the event in the target window
%   nJe      - the last index of the event in the target window
%
%   T0 ------- fTs ---------------- fTe ----- T1
%                target time window
%   1  ------- nJs ---------------- nJe ----- nNumEq
%
%   Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 1 November 2024
%   ...
%   version: 1.0.0, 1 November 2024
%
    ll_logsum = 0.0;
    for j = nJs:nJe
        s = vPar(1) + vPar(2).*sum(exp(-vPar(3).*(vCat(j,1) - vCat(1:(j-1),1)))); % \mu + A*sum
        ll_logsum = ll_logsum + log(s);
    end
    fLle = -product_hawkes(vCat,fTs,fTe,nJs,vPar) + ll_logsum; % 
end
