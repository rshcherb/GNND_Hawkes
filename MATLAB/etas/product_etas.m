function prod = product_etas(vCat,fTs,fT,nJs,vPar)
%
%   The productivity of the ETAS process without the term \mu*(fT - fTs)
%   usinng the normalization by Harte, 2010
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 23 June, 2018
%   ...
%   version 2.0.0, 2 October 2024
%
    nJ = find(vCat(:,1) <= fT,1,'last'); % find the last index for which vCat(:,1) <= fT
    vexpat = exp(vPar(4).*vCat(1:nJ,2));
    vtc = vCat(1:nJ,1)./vPar(2);
    fTsc = fTs/vPar(2);                  % fTs/c
    fTc = fT/vPar(2);                    % fTe/c
    if vPar(3) == 1.0
        prod = sum(vexpat(1:nJs-1).*log( (fTc - vtc(1:nJs-1) + 1)./(fTsc - vtc(1:nJs-1) + 1) ));
        prod = vPar(1).*vPar(2).*(prod + sum(vexpat(nJs:nJ).*log( fTc - vtc(nJs:nJ) + 1 )));
    else
        p1p = 1.0 - vPar(3);
        prod = sum(vexpat(1:nJs-1).*((fTc - vtc(1:nJs-1) + 1).^p1p - (fTsc - vtc(1:nJs-1) + 1).^p1p));
        prod = vPar(1).*vPar(2)./p1p.*(prod + sum(vexpat(nJs:nJ).*((fTc - vtc(nJs:nJ) + 1).^p1p - 1)));
    end
end
