function vMoment = mag2moment(vMag)
%
%   Function to convert magnitudes to moments
%
%   Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 17 June 2022
%   ...
%   version: 1.0.0, 17 June 2022
%
    vMoment = 10.^(1.5.*(vMag + 10.73));
end
