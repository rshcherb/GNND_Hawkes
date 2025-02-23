function vMag = moment2mag(vMoment)
%
%   Function to convert moments to magnitudes
%
%   Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 17 June 2022
%   ...
%   version: 1.0.0, 17 June 2022
%
     vMag = 2.0/3.0*log10(vMoment) - 10.73;
end
