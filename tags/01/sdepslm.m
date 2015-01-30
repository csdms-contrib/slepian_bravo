function e2lm=sdepslm(TH,L,bas)
% e2lm=SDEPSLM(TH,L,bas)
% 
% Mean-square error on the coefficients of a recovered field from
% incomplete observations (NOT observed inside a double polar cap)
% under the wrong minimization criterion, thus utter garbage, really.
%
% INPUT:
%
% TH          Angular extent of the spherical cap, in degrees
% L           Bandwidth (maximum angular degree)
% bas         Basis, 'sh' (default) or 'slepian'
% 
% OUTPUT:
%
% e2lm        Mean-square error on the coefficients
%
% What we need is the diagonal elements of the localization kernel, and
% we need them for every m so we can arrange them in pyramid form.
%
% This really is about how well a particular matrix can be represented by
% its largest eigenvectors, as it has been all along.
%
% See also SDSNEEUW, SDEPSTH, PLM2MAP
%
% Last modified by fjsimons-at-alum.mit.edu, 05/21/2009

defval('bas','sh')

% Initialize matrix
e2lm=repmat(0,L+1,L+1);
% Do this for every order
for m=0:L
  [E,V,th,C,ngl1,ngl2,K]=sdwcap2(TH,L,m,0);
  switch bas
   case 'sh'
    e2lm(m+1:L+1,m+1)=diag(K);
   case 'slepian'
    % Fill the last values in with VCUT, eps approximately
    e2lm(m+1:L+1,m+1)=[V(:) ; repmat(eps,L-m+1-length(V),1)];
  end
end

% Check the average value, which should be (1-cos(TH*pi/180))
if strcmp(bas,'sh')
  difer(abs(sum((e2lm(:,1)+2*sum(e2lm(:,2:end),2))./(2*(0:L)+1)'...
	       -(1-cos(TH*pi/180)))));
end

% Replace zeroes with NaN's for plotting
e2lm(e2lm==0)=NaN;






