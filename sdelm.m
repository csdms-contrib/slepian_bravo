function varargout=sdelm(TH,L,cp,xp)
% [lrnk,mrnk,lval,VV,Vsum]=SDELM(TH,L,cp,xp)
%
% Computes globally ranked eigenvalues for the concentration problem to
% the spherical DOUBLE POLAR CAP or the complementary LATITUDINAL BELT,
% across all angular orders m
%
% INPUT:
%
% TH      Colatitudinal radius of the DOUBLE POLAR CAP, in degrees
% L       Bandwidth
% cp      1 The double polar cap
%         2 The complementary latitudinal belt
% xp      0 Absolute orders m only
%         1 Expanded form including +/- m
% 
% OUTPUT:
% 
% lrnk    The index of the ranked eigenvalue within its own m
% mrnk    The m of the ranked eigenvalue
% lval    The sorted eigenvalues
% VV      The not globally sorted eigenvalue matrix (lrank-by-m)
% Vsum    Total sum of all the eigenvalues, including double counts
%
% Last modified by fjsimons-at-alum.mit.edu, 09.08.2005

defval('TH',40);
defval('L',18);
defval('cp',1);
defval('xp',0);

fnpl=sprintf('%s/SDELM-%i-%i-%i.mat',...
	     fullfile(getenv('IFILES'),'WIECZOREK'),...
	     TH,L,cp);

if exist(fnpl,'file')~=2
  disp('SDELM Start')
  % Initialize eigenvalue matrix
  VV=repmat(NaN,L+1,L+1);
  M=0:L;
  for m=M
    % Don't compute the eigenfunctions, just the eigenvalue
    [E,Vg,th,C,T,V]=grunbaum2(TH,L,m,0);
    switch cp
     case 1
      V=V;
     case 2
      V=1-V;
     otherwise
      error('Specify valid option: [1] Caps or [2] Belt')
    end
    % All the eigenvalues; sometimes slightly negative
    % This is only filled for possible values from m to L
    VV(1:min(L-m+1,length(V)),m+1)=V(:);
  end
  disp('SDELM End')
  
  % Figure out rank ordering for the eigenvalues, from high to low
  [a,b]=sort(VV(:));
  mrnk=repmat(0:L+1,L+1,1);
  lrnk=repmat([1:L+1]',1,L+1);
  b=b(~isnan(a)); b=flipud(b);
  a=a(~isnan(a)); a=flipud(a);
  % The ranked eigenvalues belong to this m
  mrnk=mrnk(b);
  % And they represent this number of nth eigenvalue
  lrnk=lrnk(b);
  lval=a;
  
  % Figure out total sum of the eigenvalues including double counts
  Vsum=VV;
  Vsum(:,2:end)=Vsum(:,2:end)*2;
  Vsum=sum(Vsum(~isnan(Vsum)));
  eval(sprintf('save %s lrnk mrnk lval VV Vsum',fnpl))
else
  eval(sprintf('load %s',fnpl))
end

if xp==1
  % Nearly double the amount of requested tapers for repeated abs(m)
  lval=gamini(lval,(mrnk~=0)+1);
  lrnk=gamini(lrnk,(mrnk~=0)+1);
  mrnk=gamini(mrnk,(mrnk~=0)+1);
  % Remember this was different in the old version
  mrnk((diff(mrnk)==0) & (diff(lrnk)==0))= ...
      -mrnk((diff(mrnk)==0)& (diff(lrnk)==0));
  % Take a look at the ranked orders and degrees: [lrnk(:) mrnk(:)]
  if length(lrnk)~=(L+1)^2
    error('Something wrong')
  end
  lval=lval(:);
  lrnk=lrnk(:);
  mrnk=mrnk(:);
end

% Output
varn={'lrnk','mrnk','lval','VV','Vsum'};
for index=1:nargout
  varargout{index}=eval(varn{index});
end
