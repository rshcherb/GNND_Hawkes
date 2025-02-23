function vPdf_mag = pdf_moment2mag(vMag,vPdf_mom)
%
%   Function to convert the pdf from moment domain to magnitude domain
%
%   Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 28 August 2024
%   ...
%   version: 1.0.0, 28 August 2024
%
     vPdf_mag = 1.5*log(10)*vPdf_mom.*10.^(1.5*(vMag + 10.73));
end
