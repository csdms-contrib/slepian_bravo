function lmcosip=plm2mag(lmcosi,r,a,wat,schm)
% lmcosip=PLM2MAG(lmcosi,r,a,wat,schm)
%
% Converts a (Schmidt-normalized) real-spherical-harmonic array of
% magnetic-potential cosine and sine coefficients into a same-size
% spherical harmonic array with modified coefficients, depending on given
% input and requested output, e.g. for upward and downward continuation,
% differentiation to individual components, and the like. 
%
% INPUT:
%
% lmcosi        [l m Ccos Csin] degrees, order, cos/sin potential coefficients
% r             Requested radius for the input, in m [default: a]
% a             Reference radius for the output [default: 6371200]
% wat           1 output is radial-component field [default]
%               2 output is colatitudinal tangential-component field [not done yet]
%               3 output is longitudinal tangential-component field [not done yet]
% schm          1 change normalization from being appropriate for expansion
%                 into Schmidt semi-normalized harmonics to being
%                 appropriate for expansion into 4pi-normalized ones [default] 
%                 Otherwise, keep the normalization as it is.
%
% OUPUT:
%
% lmcosip       [l m Ccos Csin] degrees, order, cos/sin coefficients of output
%
% Written by by fjsimons-at-alum.mit.edu, 03/05/2009
% Last modified by Jarno Saarimaki, 07/01/2011
% Last modified by fjsimons-at-alum.mit.edu, 03/20/2020

% Deafult input values
defval('a',6371200)
defval('r',a)
defval('wat',1)
defval('schm',1)

% Get the spherical harmonic degrees
el=lmcosi(:,1);

% Precalculate a certain upward-continuation factor
arl2=(a/r).^(el+2);

if schm==1
  % Convert coefficients from being expanded into Schmidt semi-normalized
  % to being expanded into 4pi-normalized real spherical harmonics
  lmcosi(:,3:4)=lmcosi(:,3:4)./sqrt(2*[el el]+1);
end

% Perform the upward or downward continuation portion first
lmcosi(:,3:4)=lmcosi(:,3:4).*(arl2*[1 1]);

% Convert to other things
switch wat
 case 1
  % [Blakely 1995 p. 169 after eq. (8.20)]
  disp('Calculating radial field from potential')
  lmcosip(:,3:4)=lmcosi(:,3:4).*((el+1)*[1 1]);
 case 2
  % Kono 
  disp('Calculating colatitudinal tangential field from potential')
 case 3
  % Kono
  disp('Calculating longitudinal tangential field from potential')
end
