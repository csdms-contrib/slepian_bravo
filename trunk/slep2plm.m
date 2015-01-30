function varargout=slep2plm(falpha,TH,L,phi,theta,omega)
% [lmcosi,V,N]=SLEP2PLM(falpha,TH,L,phi,theta,omega)
%
% Finds the expansion coefficiens into a SINGLE-CAP, potentially rotated,
% Slepian basis of a function whose real spherical harmonic expansion
% coefficients are known.
%
% INPUT:
%
% falpha     Slepian expansion coefficients
% TH         Radius of the concentration region (degrees) OR
%            'england', 'eurasia',  'namerica', 'australia', 'greenland'
%            'africa', 'samerica', 'amazon', 'orinoco'
% L          Bandwidth of the window [default: bandwidth of the data]
% phi        Longitude of the center of the tapers (degrees)
% theta      Colatitude of the center of the tapers (degrees)
% omega      Anticlockwise azimuthal rotation of the tapers (degrees)
%
% OUTPUT:
%
% lmcosi     Standard-type real spherical harmonic expansion coefficients
% V          The eigenvalues of the Slepian functions in question
% N          The Shannon number
%
% See also: PTOSLEP, GLMALPHA, GLMALPHAPTO
%
% Last modified by fjsimons-at-alum.mit.edu, 06/27/2012

% Supply defaults
defval('TH',30)
defval('L',18)
defval('phi',0)
defval('theta',0)
defval('omega',0)
defval('nosort',0)

% If it is the standard North-Polar cap, it's easy
if phi==0 && theta==0 && omega==0
  % Get the Slepian basis; definitely not block-sorted as for the rotated
  % versions this will make no sense at all anymore
  [G,V,EL,EM,N,GM2AL,MTAP,IMTAP]=glmalpha(TH,L,[],0);
else
  % Need to get a complete GLMALPHA but for the rotated basis
  % Definitely, "single-order" has lost its meaning here, but the MTAP
  % will still identify what the order of the unrotated original was
  [G,V,EL,EM,N,GM2AL,MTAP,IMTAP]=glmalphapto(TH,L,phi,theta,omega);
end

% Sort by decreasing eigenvalue, rely on the falphas to be similarly sorted
[V,vi]=sort(V,'descend');
G=G(:,vi); if ~isnan(MTAP); MTAP=MTAP(vi); end

% Now perform the inverse expansion
% Get the mapping from LMCOSI into not-block-sorted GLMALPHA
[~,~,~,lmcosi,~,~,~,~,~,ronm]=addmon(L);

% Perform the expansion of the signal into the Slepian basis
% and stick these coefficients in at the right places
lmcosi(2*size(lmcosi,1)+ronm(1:(L+1)^2))=G*falpha;

% Collect output
varns={lmcosi,V,N,MTAP};
varargout=varns(1:nargout);
