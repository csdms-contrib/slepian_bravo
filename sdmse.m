function [eta,mseO,mseR,mseC,fO,fR,fC]=sdmse(SN,eta,TH,L,phi)
% [eta,mseO,mseR,mseC,fO,fR,fC]=sdmse(SN,eta,TH,L,phi)
%
% Calculates the colatitudinal mean square error AVERAGED over the 
% (complementary) region of observation or over the entire sphere
% for spherical harmonic expansions, by Gauss-Legendre averaging
%
% INPUT:
%
% SN       Signal-to-noise ratio, scalar
% eta      Damping parameter (0 <= eta <= 1) or 
%          A single number identifiying a linearly spaced set
% TH       Angular extent of the spherical caps, in degrees
% L        Bandwidth (maximum angular degree)
% phi      Perhaps you want to put in a particular longitude (default: not)
%
% OUTPUT:
%
% eta      Damping parameter
% mseO     Mean square error averaged over the sphere
% mseR     Mean square error averaged over the region
% mseC     Mean square error averaged over the complementary region 
% fO       Eta of minimum mean square error over the sphere 
% fR       Eta of minimum mean square error over the region 
% fC       Eta of minimum mean square error over the complementary region 
%
% Last modified by fjsimons-at-alum.mit.edu, 09/23/2023
% 
% See also SDETA, SDMSK, SDMSE2

defval('SN',10)
defval('eta',100)
defval('TH',10)
defval('L',45)
defval('phi',NaN);

% Signal and noise variances
S=1;
N=S/SN;

% (Number of) Damping or truncation levels tried
if eta>1 & length(eta)==1
  eta=linspace(0,1,eta);
end

% For the double polar cap
sord=2;
% Start doing colatitudes only
% Gauss-Legendre points ! It's a sum of squared polynomials...
intvO=[-1 1];
[wO,xO]=gausslegendrecof(2*L,[],[-1 1]);
intvR=cos([180-TH TH]*pi/180);
[wR,xR]=gausslegendrecof(2*L,[],intvR);
intvC=[cos(TH*pi/180) 1];
[wC,xC]=gausslegendrecof(2*L,[],intvC);
% Specify integration points
thetaO=acos(xO);
thetaR=acos(xR);
thetaC=acos(xC);

% The area of measurement is the belt
srt='belt';

% Load precomputed Slepian functions?
fnpl=sprintf('%s/SDERRGL-%i-%i-%i.mat',...
	     fullfile(getenv('IFILES'),'SDERR'),TH,L,phi);

if exist(fnpl,'file')==2 
  eval(sprintf('load %s',fnpl))
  disp(sprintf('load %s',fnpl))
else
  % Calculate eigenfunctions and eigenvalues
  [GO,Vo,jk1,jk2,jk3,KA,K]=galpha(TH,L,sord,thetaO,phi,srt);
  % Calculate eigenfunctions and eigenvalues
  GR=galpha(TH,L,sord,thetaR,phi,srt);
  % Calculate eigenfunctions and eigenvalues
  GC=galpha(TH,L,sord,thetaC,phi,srt);

  eval(sprintf(['save %s '...
		'GO GR GC Vo thetaO thetaR thetaC KA K'],...
		fnpl))
end

% Load precomputed variance results?
fnpl=sprintf('%s/SDERRGL-%i-%i-%i-%i-%i.mat',...
	     fullfile(getenv('IFILES'),'SDERR'),SN,length(eta),TH,L,phi);

if exist(fnpl,'file')==2 & eta(1)==0 & eta(end)==1 
  eval(sprintf('load %s',fnpl))
  disp(sprintf('load %s',fnpl))
else
  % DAMPED SPHERICAL HARMONIC INVERSIONS
  mseO=repmat(NaN,1,length(eta));
  mseR=repmat(NaN,1,length(eta));
  mseC=repmat(NaN,1,length(eta));
  
  for index=1:length(eta)
    % Variance matrix
    V=sdvar(N,eta(index),Vo);
    % Bias matrix
    B=sdbias(S,eta(index),Vo);
    % Error matrix for the damped spherical harmonic case
    % Diagonal terms only, at identical locations
    errSHO=diag(GO'*(V+eta(index)^2*B)*GO);
    errSHR=diag(GR'*(V+eta(index)^2*B)*GR);
    errSHC=diag(GC'*(V+eta(index)^2*B)*GC);

    % Mean square error over the entire sphere
    % Integral over the entire sphere, divided by area
    mseO(index)=wO(:)'*errSHO/range(intvO);
    % Integral over the belt
    mseR(index)=wR(:)'*errSHR/range(intvR);
    % Integral over both caps (factor two cancels)
    mseC(index)=wC(:)'*errSHC/range(intvC);
  end

  % Check that these things add up correctly
  difer(mseO-mseR*range(intvR)/2-2*mseC*range(intvC)/2,[],0);

  % So much for the calculation, now for the prediction
  for und=1:length(eta)
    da=Vo+eta(und)*(1-Vo);
    fO(und)=eta(und)-sum(Vo.*(1-Vo)./da.^3)...
	   ./sum(Vo.*(1-Vo).^2./da.^3)/SN;
    fR(und)=eta(und)-sum(Vo.^2.*(1-Vo)./da.^3)...
	    ./sum(Vo.^2.*(1-Vo).^2./da.^3)/SN;
    fC(und)=eta(und)-sum(Vo.*(1-Vo).^2./da.^3)...
	    ./sum(Vo.*(1-Vo).^3./da.^3)/SN;
  end

  % Find minimum variance and compare to theory
  [a,fO]=min(abs(fO));
  difer(mseO(fO)-min(mseO))
  fO=eta(fO);
  [c,fR]=min(abs(fR)); 
  difer(mseR(fR)-min(mseR))
  fR=eta(fR);
  [e,fC]=min(abs(fC)); 
  difer(mseC(fC)-min(mseC))
  fC=eta(fC);
  
  if eta(1)==0 & eta(end)==1
    eval(sprintf([...
	'save %s eta mseO mseR mseC fO fR fC'],fnpl))
  end
end




