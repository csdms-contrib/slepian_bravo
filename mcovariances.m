function [Sllp,l]=mcovariances(L,TH,sord,Lmax)
% [Sllp,l]=MCOVARIANCES(L,TH,sord,Lmax)
%
% Calculates the eigenvalue-weighted multitaper covariance matrix
%
% INPUT:
%
% L       Bandwidth of the multitaper windows, always a scalar
% TH      Colatitudinal radius of the cap, in degrees <=180, may be vector
%         Colatitudinal halfwidth of the cut, degrees <=90, may be vector
% sord    1 Single cap of diameter 2TH [default], always a scalar
%         2 Double cap left by subtracting belt of width 2TH, always scalar
% Lmax    Maximum spherical harmonic degree
%
% OUTPUT:
%
% Sllp    The covariance matrix indexed by degree l and l' 
% l       The degrees at which this is evaluated [0->Lmax]
%
% Last modified by fjsimons-at-alum.mit.edu, 03/07/2007

% Plot with imagefnan

defval('L',10)
defval('TH',80)
defval('sord',2)
defval('Lmax',100)
defval('xver',1)

% A quick routine, first in line for improvement

% Load the Gammap matrix
[Gp,p]=gammap(L,TH,sord,0);

% Load the ZEROJ database once (don't ask for a particular database by
% bandwidth, but let Matlab determine the best available database
[jk,C0,S0,Leff]=zeroj(Lmax+2*L,0,0);

% Initialize
Sllp=nan(Lmax+1,Lmax+1);

h=waitbar(0,'Construction covariance matrix');
for l=0:Lmax
  for lp=0:l
    Sllp(l+1,lp+1)=sum((2*p+1).*Gp.*...
	    zeroj(repmat(l,1,length(p)),p,repmat(lp,1,length(p)),...
		  Leff,[],C0,S0).^2);
    Sllp(lp+1,l+1)=Sllp(l+1,lp+1);
    waitbar(l/Lmax,h)
  end
end
Sllp=Sllp/2/pi;
close(h)

% Get the degrees back
l=0:Lmax;

if xver==1
  % For comparison with mvarratio transform to the ratio to the WS
  difer(mvarratios(L,TH,sord,0:Lmax,0)'-diag(Sllp).*(2*l(:)+1)/2)
end

keyboard
