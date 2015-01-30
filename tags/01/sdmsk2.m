function [kay,mskO,mskR,mskC,kO,kR,kC,K,varO,bs2O]=sdmsk2(SN,kay,TH,L)
% [kay,mskO,mskR,mskC,kO,kR,kC,K,varO,bs2O]=sdmsk2(SN,kay,TH,L)
%
% Calculates the mean square error AVERAGED over the 
% (complementary) region of observation or over the entire sphere
% for truncated Slepian expansions, by hand.
%
% INPUT:
%
% SN          Signal-to-noise ratio, scalar
% kay         Truncation parameter (0 <= kay <= (L+1)^2) (default: all)
% TH          Angular extent of the spherical caps, in degrees
% L           Bandwidth (maximum angular degree)
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
% varO        The variance component of mskO 
% bs2O        The bias^2 component of mskO 
%
% Last modified by fjsimons-at-alum.mit.edu, 04/13/2007
% 
% See also SDKAY, SDMSE, SDMSK

defval('SN',10)
defval('kay',[0:(L+1)^2])
defval('TH',10)
defval('L',45)
defval('ntheta',720)
theta=linspace(0,pi,720);
defval('phi',NaN);

% Signal and noise variances
S=1;
N=S/SN;

% For the double polar cap
sord=2;
% The area of measurement is the belt
srt='belt';

% Area element for single or double cap
A=2*pi*(1-cos(TH*pi/180));
if sord==2
  A=2*A;
end
 
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

% Load precomputed variance results?
fnpl=sprintf('%s/SDERRGLK2-%i-%i-%i-%i.mat',...
	     fullfile(getenv('IFILES'),'SDERR'),SN,(L+1)^2+1,TH,L);

if exist(fnpl,'file')==2 & kay(1)==0 & kay(end)==(L+1)^2
  eval(sprintf('load %s',fnpl))
  disp(sprintf('load %s',fnpl))
else
  % TRUNCATED SLEPIAN APPROACH
  mskO=repmat(NaN,1,length(kay));
  mskR=repmat(NaN,1,length(kay));
  mskC=repmat(NaN,1,length(kay));
  varO=repmat(NaN,1,length(kay));
  bs2O=repmat(NaN,1,length(kay));
  % Bias and variance terms
  V=1./Vo;
  B=repmat(1,size(Vo));
  for index=1:length(kay)
    % kay-indices
    leK=1:kay(index);
    gtK=kay(index)+1:(L+1)^2;
    % Bias and variance for real
    varO(index)=N*sum(V(leK))/4/pi;
    bs2O(index)=S*sum(B(gtK))/4/pi;
    % This is the mean square error averaged over the area...
    % ... of the sphere
    mskO(index)=varO(index)+bs2O(index);
    % ... of the belt
    mskR(index)=(N*sum(B(leK))+S*sum(B(gtK).*Vo(gtK)))/(4*pi-A);
    % ... of the polar caps
    mskC(index)=(N*sum(V(leK).*(1-Vo(leK)))+S*sum(B(gtK).*(1-Vo(gtK))))/A;
  end

  % Check that these things add up correctly
  difer(mskO-mskR*(1-A/4/pi)-mskC*A/4/pi,[],0);

  % So much for the calculation, now for the prediction
  [jk,kO]=min(abs(Vo-1/SN));
  kR=kO; kC=kO;
  % Find minimum variance and compare to theory; note kay(kO+1)=kO
  difer(mskO(kO+1)-min(mskO),[],0)

  % Use the real minimum of the variance, shall we?
  [jk,kO]=min(mskO);
  kO=kay(kO);
  difer(mskR(kR)-min(mskR),[],0)
  [jk,kR]=min(mskR);
  kR=kay(kR);
  difer(mskC(kC)-min(mskC),[],0)
  [jk,kC]=min(mskC);
  kC=kay(kC);
  
  if kay(1)==0 & kay(end)==(L+1)^2
    eval(sprintf([...
	'save %s kay mskO mskR mskC kO kR kC K varO bs2O'],fnpl))
  end
end




