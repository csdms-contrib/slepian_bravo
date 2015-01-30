function [mskTH,theta,varkTH,bsk2TH]=sderk(SN,TH,L,kay,ntheta,phi)
% [mskTH,theta,varkTH,bsk2TH]=SDERK(SN,TH,L,kay,ntheta,phi)
%
% Calculates the colatitudinal mean square error, NOT AVERAGED
% for a particular Slepian truncation level.
%
% INPUT:
%
% SN          Signal-to-noise ratio
% TH          Angular extent of the spherical caps, in degrees
% L           Bandwidth (maximum angular degree)
% kay         The truncation factor
% ntheta      Number of colatitudes (default: 720)
% phi         Perhaps you want to put in a particular longitude (default: not)
%
% OUTPUT:
%
% mskTH       Mean square error at requested colatitudes
% theta       Colatitudes used, linearly spaced values
% varkTH       Variance contributiuon to the mse
% bs2kTH       Bias contributiuon to the mse
% 
% Last modified by fjsimons-at-alum.mit.edu, 04/13/2007
% 
% See also SDMSK

defval('SN',10)
defval('TH',10)
defval('L',45)
defval('ntheta',720)
defval('phi',NaN);

% Signal and noise variances
S=1;
N=S/SN;

% For the double polar cap
sord=2;

% Start doing colatitudes only
theta=linspace(0,pi,720);

% The area of measurement is the belt
srt='belt';

% Precompute the Slepian solutions
fnpl=sprintf('%s/SDERR-%i-%i-%i-%i.mat',...
	     fullfile(getenv('IFILES'),'SDERR'),TH,L,ntheta,phi);

if exist(fnpl,'file')==2
  eval(sprintf('load %s',fnpl))
  disp(sprintf('load %s',fnpl))
else
  % Calculate eigenfunctions and eigenvalues
  [G,Vo,jk1,jk2,jk3,KA,K]=galpha(TH,L,sord,theta,phi,srt);
  eval(sprintf(['save %s G Vo theta K'],fnpl))
end

fnpl=sprintf('%s/SDERRSG-%i-%i-%i-%i-%i.mat',...
	     fullfile(getenv('IFILES'),'SDERR'),SN,TH,L,kay,phi);

% Excessive verification?
xver=1;
if exist(fnpl,'file')==2
  eval(sprintf('load %s',fnpl))
  disp(sprintf('load %s',fnpl))
else
  % TRUNCATED SPHERICAL HARMONIC INVERSIONS
  leK=1:kay;
  gtK=kay+1:(L+1)^2;
  % Variance matrix
  V=sdvar(N,NaN,Vo(leK));
  % Bias matrix
  B=sdbias(S,NaN,Vo(gtK));
  % Error matrix for the truncated Slepian case
  % Diagonal terms only, at identical locations
  % mskTH=diag(G(leK,:)'*V*G(leK,:)+...
  %	     G(gtK,:)'*B*G(gtK,:));
  varkTH=diag(G(leK,:)'*V*G(leK,:));
  bsk2TH=diag(G(gtK,:)'*B*G(gtK,:));
  mskTH=varkTH+bsk2TH;
  % Compare this! Pointless, though
  if xver==1
    V=1./Vo;
    B=repmat(1,size(Vo));
    varTH=[N*V(leK)*(G(leK,:).^2)]';
    bs2TH=[S*B(gtK)*(G(gtK,:).^2)]';
    msTH=varTH+bs2TH;
    difer(msTH-mskTH)
    difer(bs2TH-bsk2TH)
    difer(varTH-varkTH)
  end
  eval(sprintf(['save %s mskTH theta varkTH bsk2TH'],fnpl))
end





