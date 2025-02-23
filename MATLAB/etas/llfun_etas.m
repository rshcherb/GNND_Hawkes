function fLle = llfun_etas(vPar,vCat,fTs,fTe,nJs,nJe)
%
%   The log-likelihood function for the ETAS conditional rate
%   using the normalization by Harte (2010)
%
%   vPar     - the parameters of the ETAS model as variables: \mu = vPar(1); K = vPar(2); c = vPar(3); p = vPar(4); \alpha = vPar(5);
%   vCat      - earthquake times and magnitudes with the subtracted reference magnitude m_0 (vTM(:,2) - m_0)
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
%   version: 1.0.0, February 19, 2015
%   ...
%   version: 4.1.1, 20 November 2024
%
    vexpat = exp(vPar(5).*vCat(1:nJe,2));
    vtc = vCat(1:nJe,1)./vPar(3);

    ll_logsum = 0.0;
    for j = nJs:nJe
        s = vPar(1) + vPar(2).*sum(vexpat(1:(j-1))./((vtc(j) - vtc(1:(j-1)) + 1).^vPar(4))); % \mu + K*sum
        ll_logsum = ll_logsum + log(s);
    end
    fLle = -vPar(1).*(fTe-fTs) - product_etas(vCat,fTs,fTe,nJs,vPar(2:5)) + ll_logsum; % 
end
