function [mseTH,theta,varTH,bs2TH]=sderr(SN,TH,L,eta,ntheta,phi)
% [mseTH,theta,varTH,bs2TH]=SDERR(SN,TH,L,eta,ntheta,phi)
%
% Calculates the colatitudinal, not averaged, mean square estimation
% error of the geodetic estimation problem at altitude a=0, for a
% particular polar gap, truncation bandwidth, and a single damping level
% used in the damped spherical harmonic solution approach.
%
% INPUT:
%
% SN          Signal-to-noise ratio (default: 10)
% TH          Angular extent of both spherical caps, in degrees (default: 10)
% L           Bandwidth, maximum angular degree (default: 45)
% eta         The damping factor
% ntheta      Number of colatitudes (default: 720)
% phi         Perhaps you want to put in a particular longitude (default: not)
%
% OUTPUT:
%
% mseTH       Mean square error at requested colatitudes
% theta       Colatitudes used, linearly spaced values
% varTH       Variance contributiuon to the mse
% bs2TH       Bias contributiuon to the mse
%
% Last modified by fjsimons-at-alum.mit.edu, 04/13/2007
% 
% See also SDMSE, SDERRS, SDVAR, SDBIAS

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
theta=linspace(0,pi,ntheta);

% The area of measurement is the belt
srt='belt';

% Precompute the Slepian solutions
fnpl=sprintf('%s/SDERR-%i-%i-%i-%i.mat',...
	     fullfile(getenv('IFILES'),'SDERR'),TH,L,ntheta,phi);

if exist(fnpl,'file')==2 
  eval(sprintf('load %s',fnpl))
else
  % Calculate eigenfunctions and eigenvalues
  [G,Vo,jk1,jk2,jk3,KA,K]=galpha(TH,L,sord,theta,phi,srt);
  eval(sprintf(['save %s G Vo theta K'],fnpl))
end

fnpl=sprintf('%s/SDERRSH-%i-%i-%i-%5.3f-%i.mat',...
	     fullfile(getenv('IFILES'),'SDERR'),SN,TH,L,eta,phi);

if exist(fnpl,'file')==2
  eval(sprintf('load %s',fnpl))
  disp(sprintf('load %s',fnpl))
else
  % DAMPED SPHERICAL HARMONIC INVERSIONS
  mseTH=repmat(NaN,length(theta),length(eta));
  
  % Variance matrix
  V=sdvar(N,eta,Vo);
  % Bias matrix
  B=sdbias(S,eta,Vo);
  % Error matrix for the damped spherical harmonic case
  % Diagonal terms only, at identical locations
  % mseTH=diag(G'*(V+eta^2*B)*G);
  varTH=diag(G'*V*G);
  % Speed this up if eta=0
  if eta==0
    bs2TH=0;
  else
    bs2TH=eta^2*diag(G'*B*G);
  end
  mseTH=varTH+bs2TH;
  eval(sprintf(['save %s mseTH theta varTH bs2TH'],fnpl))
end





