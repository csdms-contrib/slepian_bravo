function [kay,mskO,mskR,mskC,kO,kR,kC,K]=sdmsk(SN,kay,TH,L,phi)
% [kay,mskO,mskR,mskC,kO,kR,kC,K]=sdmsk(SN,kay,TH,L,phi)
%
% Calculates the mean square error AVERAGED over the 
% (complementary) region of observation or over the entire sphere
% for truncated Slepian expansions, by Gauss-Legendre averaging
%
% INPUT:
%
% SN          Signal-to-noise ratio, scalar
% kay         Truncation parameter (0 <= kay <= (L+1)^2) (default: all)
% TH          Angular extent of the spherical caps, in degrees
% L           Bandwidth (maximum angular degree)
% phi         Perhaps you want to put in a particular longitude (default: not)
%
% OUTPUT:
%
% kay         Truncation levels
% mskO        Mean square error averaged over the sphere
% mskR        Mean square error averaged over the region
% mskC        Mean square error averaged over the complementary region 
% kO          Kay of minimum mean square error over the sphere 
% kR          Kay of minimum mean square error over the region 
% kC          Kay of minimum mean square error over the complementary region 
% K           Shannon number of this problem
%
% Last modified by fjsimons-at-alum.mit.edu, 04/13/2007
% 
% See also SDKAY, SDMSE

defval('SN',10)
defval('kay',[0:(L+1)^2])
defval('TH',10)
defval('L',45)
defval('phi',NaN);

% Signal and noise variances
S=1;
N=S/SN;

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
% Longitude: if NaN do not include the cos(m phi) or sin(m phi)

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
  [GO,Vo,EM,jk2,jk3,KA,K]=galpha(TH,L,sord,thetaO,phi,srt);
  % Calculate eigenfunctions and eigenvalues
  GR=galpha(TH,L,sord,thetaR,phi,srt);
  % Calculate eigenfunctions and eigenvalues
  GC=galpha(TH,L,sord,thetaC,phi,srt);

  eval(sprintf(['save %s '...
		'GO GR GC Vo thetaO thetaR thetaC KA K'],...
		fnpl))
end

% Load precomputed variance results?
fnpl=sprintf('%s/SDERRGLK-%i-%i-%i-%i-%i.mat',...
	     fullfile(getenv('IFILES'),'SDERR'),SN,(L+1)^2+1,TH,L,phi);

if exist(fnpl,'file')==2 & kay(1)==0 & kay(end)==(L+1)^2
  eval(sprintf('load %s',fnpl))
  disp(sprintf('load %s',fnpl))
else
  % TRUNCATED SLEPIAN APPROACH
  mskO=repmat(NaN,1,length(kay));
  mskR=repmat(NaN,1,length(kay));
  mskC=repmat(NaN,1,length(kay));
  
  for index=1:length(kay)
    leK=1:kay(index);
    gtK=kay(index)+1:(L+1)^2;
    % Variance matrix
    V=sdvar(N,NaN,Vo(leK));
    % Bias matrix
    B=sdbias(S,NaN,Vo(gtK));
    % Error matrix for the truncated Slepian case
    % Diagonal terms only, at identical locations
    errSGO=diag(GO(leK,:)'*V*GO(leK,:)+...
		GO(gtK,:)'*B*GO(gtK,:));
    errSGR=diag(GR(leK,:)'*V*GR(leK,:)+...
		GR(gtK,:)'*B*GR(gtK,:));
    errSGC=diag(GC(leK,:)'*V*GC(leK,:)+...
		GC(gtK,:)'*B*GC(gtK,:));

    % Mean square error over the entire colatitudinal range
    % Integral over the entire sphere, divided by area
    mskO(index)=wO(:)'*errSGO/range(intvO);
    % Integral over the belt
    mskR(index)=wR(:)'*errSGR/range(intvR);
    % Integral over both caps (factor two cancels)
    mskC(index)=wC(:)'*errSGC/range(intvC);
  end

  % Check that these things add up correctly
  difer(mskO-mskR*range(intvR)/2-2*mskC*range(intvC)/2,[],0);

  % So much for the calculation, now for the prediction
  [jk,kO]=min(abs(Vo-1/SN));
  kR=kO;
  kC=kO;

  % Find minimum variance and compare to theory
  difer(mskO(kO)-min(mskO),[],0)
  kO=kay(kO);
  difer(mskR(kR)-min(mskR),[],0)
  kR=kay(kR);
  difer(mskC(kC)-min(mskC),[],0)
  kC=kay(kC);
  
  if kay(1)==0 & kay(end)==(L+1)^2
    eval(sprintf([...
	'save %s kay mskO mskR mskC kO kR kC K'],fnpl))
  end
end
