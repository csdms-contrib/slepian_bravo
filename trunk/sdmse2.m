function [eta,mseO,mseR,mseC,fO,fR,fC,varO,bs2O]=sdmse2(SN,eta,TH,L)
% [eta,mseO,mseR,mseC,fO,fR,fC,varO,bs2O]=sdmse2(SN,eta,TH,L)
%
% Calculates the colatitudinal mean square error AVERAGED over the 
% (complementary) region of observation or over the entire sphere,
% calculated analytically, not by integration.
%
% INPUT:
%
% SN         Signal-to-noise ratio, scalar
% eta        Damping parameter (0 <= eta <= 1) or 
%            A single number identifiying a linearly spaced set
% TH         Angular extent of the spherical caps, in degrees
% L          Bandwidth (maximum angular degree)
%
% OUTPUT:
%
% eta        Damping parameter
% mseO       Mean square error averaged over the sphere
% mseR       Mean square error averaged over the region
% mseC       Mean square error averaged over the complementary region 
% fO         Eta of minimum mean square error over the sphere 
% fR         Eta of minimum mean square error over the region 
% fC         Eta of minimum mean square error over the complementary region 
% varO       The variance component of mseO 
% bs2O       The bias^2 component of mseO 
%
% Last modified by fjsimons-at-alum.mit.edu, 04/13/2007
% 
% See also SDETA, SDMSK, SDMSE

defval('SN',10)
defval('eta',100)
defval('TH',10)
defval('L',45)
defval('ntheta',720)
defval('phi',NaN);

theta=linspace(0,pi,720);
% Signal and noise variances
S=1;
N=S/SN;

% (Number of) Damping or truncation levels tried
if eta>1 & length(eta)==1
  eta=linspace(0,1,eta);
end

% For the double polar cap
sord=2;
% The area of measurement is the belt
srt='belt';

% Area element
A=2*pi*(1-cos(TH*pi/180));
if sord==2
  A=2*A;
end

% Load precomputed Slepian functions?
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

% Load precomputed variance results?
fnpl=sprintf('%s/SDERRGL2-%i-%i-%i-%i.mat',...
	     fullfile(getenv('IFILES'),'SDERR'),SN,length(eta),TH,L);

if exist(fnpl,'file')==2 & eta(1)==0 & eta(end)==1
  eval(sprintf('load %s',fnpl))
  disp(sprintf('load %s',fnpl))
else
  % DAMPED SPHERICAL HARMONIC INVERSIONS
  mseO=repmat(NaN,1,length(eta));
  mseR=repmat(NaN,1,length(eta));
  mseC=repmat(NaN,1,length(eta));
  varO=repmat(NaN,1,length(eta));
  bs2O=repmat(NaN,1,length(eta));
  
  for index=1:length(eta)
    % Bias and variance terms
    V=Vo./[Vo+eta(index)*(1-Vo)].^2;
    B=(1-Vo).^2./[Vo+eta(index)*(1-Vo)].^2;
    % Bias and variance for real
    varO(index)=N*sum(V)/4/pi;
    bs2O(index)=S*eta(index)^2*sum(B)/4/pi;
    % This is the mean square error averaged over the area...
    % ... of the sphere
    mseO(index)=varO(index)+bs2O(index);
    % ... of the belt
    mseR(index)=(N*sum(V.*Vo)+S*eta(index)^2*sum(B.*Vo))/(4*pi-A);
    % ... of the polar caps
    mseC(index)=(N*sum(V.*(1-Vo))+S*eta(index)^2*sum(B.*(1-Vo)))/A;
  end

  % Check that these things add up correctly
  difer(mseO-mseR*(1-A/4/pi)-mseC*A/4/pi,[],0);

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
	'save %s eta mseO mseR mseC fO fR fC varO bs2O'],fnpl))
  end
end




