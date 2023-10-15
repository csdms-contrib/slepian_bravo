function [Sllp,l]=mcovariances(Sl,l,L,TH,sord,Lmax)
% [Sllp,l]=MCOVARIANCES(Sl,l,L,TH,sord,Lmax)
%
% After calculating a multitaper power spectral density, calculates the
% appropriate error bars based on it, according to Dahlen & Simons
% (2008). Not for whole-sphere data, for which the simple formula in 
% Dahlen & Simons eq. (47) applies, see also RB VII p. 32.
%
% INPUT:
%
% Sl      The power spectral density whose variance we are calculating
% l       Degrees at which the calculation is to be carried out 
% L       Bandwidth of the tapers considered
% TH      Colatitudinal radius of the cap, in degrees <=180, may be vector
%         Colatitudinal halfwidth of the cut, degrees <=90, may be vector
% sord    1 Single cap of diameter 2TH [default], always a scalar
%         2 Double cap left by subtracting belt of width 2TH, always scalar
%
% OUTPUT:
%
% Sllp    The covariance matrix indexed by degree l and l' 
% l       The degrees at which this is evaluated [0->Lmax]
%
% Last modified by fjsimons-at-alum.mit.edu, 11/17/2011

% Plot with imagefnan

defval('L',10)
defval('TH',80)
defval('sord',2)
defval('xver',1)

% Load the Gammap matrix
[Gp,p]=gammap(L,TH,sord,0);
o
% Highest degree of a zero-j database, else, switch to approximation
zjmax=500;

% The maximum degree used from the zero-j database under this construction
Lmax=min(max(l),zjmax);

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
		  Leff,[],C0,S0).^2)*Sl(l+1)*Sl(lp+1);
    Sllp(lp+1,l+1)=Sllp(l+1,lp+1);
    waitbar(l/Lmax,h)
  end
end
Sllp=Sllp/2/pi;
close(h)

if xver==1
  % For comparison with mvarratio transform to the ratio to the WS
  difer(mvarratios(L,TH,sord,0:Lmax,0)'-diag(Sllp).*(2*l(:)+1)/2)
end

keyboard
