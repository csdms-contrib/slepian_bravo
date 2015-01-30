function sdgelderen2
% SDGELDEREN2
%
% Simons & Dahlen (2005)
%
% Comparison of SDGELDEREN with SDSUMALL2 on the same plot.
%
% Last modified by fjsimons-at-alum.mit.edu, 05/21/2009

% First, run this
sdgelderen

defval('TH',[3 5 7 10])
defval('ZL',[90 180 360]);
nth=720;

for ondex=1:length(TH)
  for index=1:length(ZL)
    L=ZL(index);
    % Calculate the Shannon number for the double cap
    N(ondex)=(L+1)^2*(1-cos(TH(ondex)*pi/180))/2*2;
    
    % Perhaps we had this saved already?
    fnpl=sprintf('%s/SDSUMALL2-%i-%i-%i.mat',...
		 fullfile(getenv('IFILES'),'WIECZOREK'),...
		 TH(ondex),L,(L+1)^2);
    
    % Calculate cumulative sums weighted by eigenvalue!
    if exist(fnpl,'file')~=2 
      F=0; G=0;
      for m=0:L
	% Just calculate at a single longitude
	[E,Vg,th,C,T,V]=grunbaum2(TH(ondex),L,m,nth,1);
	% There is no factor of two - sum over cos^2 and sin^2 is one
	% Weighted by the eigenvalue
	F=F+V*(E.^2)';
	% Not weighted by the eigenvalue
	G=G+ones(size(V))*(E.^2)';
      end
      save(fnpl,'F','G')
    end
  end
end

% Now plot them
cols=grey(5);
ah=flipud(getkids(gcf));
A4p=(L+1)^2/4/pi;
for ondex=1:length(TH)
  axes(ah(ondex))
  hold on
  fnpl=sprintf('%s/SDSUMALL2-%i-%i-%i.mat',...
	       fullfile(getenv('IFILES'),'WIECZOREK'),...
	       TH(ondex),L,(L+1)^2);
  eval(sprintf('load %s',fnpl))
  p(:,ondex)=plot(linspace(0,180,length(F)),F*100,'Color',cols);
  hold on
end

fig2print(gcf,'portrait')
figdisp





