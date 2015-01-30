function lmcosip=plm2mag(lmcosiS,r,a,wat)
% lmcosip=plm2mag(lmcosiS,r,a,wat)
%
% INPUT:
%
% lmcosiS       [l m Ccos Csin] degrees, order, coefficients in the
%               Schmidt semi-normalized basis
% r             Requested radius, in m [default: a]
% a             Reference radius [default: 6371200]
% wat           1 radial-component magnetic field u_r [default]
%               2 colatitudinal tangential-component u_th
%               3 longitudinal tangential-component u_ph
%
% OUPUT:
%
% lmcosip       [l m Ccos Csin] degrees, order, coefficients of output in
%               the 4pi fully normalized basis good for PLM2XYZ
%
% Last modified by fjsimons-at-alum.mit.edu, 07/13/2012

defval('a',6371200)
defval('r',a)
defval('wat',1)

% Get the spherical harmonic degree
el=lmcosiS(:,1);
arl2=(a/r).^(el+2);

% Convert from Schmidt semi-normalized to 4pi-normalized
lmcosi(:,3:4)=lmcosiS(:,3:4)./sqrt(2*[el el]+1);
lmcosip=lmcosiS;

% Convert
switch wat
 case 1
  % Blakely p. 169 (8.20)
  fact=(el+1).*arl2*[1 1];
  disp('Calculating radial field from potential')
 case 2
  % Kono 
  disp('Calculating colatitudinal tangential field from potential')
 case 3
  % Kono
  disp('Calculating longitudinal tangential field from potential')
end

% Factor in the factor
lmcosip(:,3:4)=fact.*lmcosiS(:,3:4);
