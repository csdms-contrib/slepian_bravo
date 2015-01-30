function sdsumall2
% SDSUMALL2
%
% Simons & Dahlen and Wieczorek (2005)
% Sums ALL eigenfunctions for small regions and large bandwidths
% of DOUBLE POLAR CAPS, weighted by the eigenvalue. Times 100.
%
% See also SDSUMALL, of which this is a simplified version.
% See also SDSUMK, tested to be equivalent.
%
% Last modified by fjsimons-at-alum.mit.edu, 08.08.2008

TH=[3 5 7 10]; 
ZL=[45 90 180];
ZL=90;

for ondex=1:length(TH)
  for index=1:length(ZL)
    L=ZL(index);
    legso{index}=sprintf('L = %i',L);
    % Why so low? This may be a problem.
    nth=L+2;
    % Calculate the Shannon number for the double cap
    N(ondex)=(L+1)^2*(1-cos(TH(ondex)*pi/180))/2*2;

    % Perhaps we had this saved already?
    fnpl=sprintf('%s/SDSUMALL2-%i-%i-%i.mat',...
		 fullfile(getenv('IFILES'),'WIECZOREK'),...
		 TH(ondex),L,(L+1)^2);
    
    % Calculate cumulative sums weighted by eigenvalue!
    if exist(fnpl,'file')~=2 
      F=0; G=0; GG=0;
      % Do every order just once since we don't do longitudes
      for m=0:L
	% Just calculate at a single longitude here
	[E,Vg,th,C,T,V]=grunbaum2(TH(ondex),L,m,nth,1);
	% Single longitude, but the sqrt(2) is already in there!
	F=F+V*(E.^2)';
	% Not weighted by the eigenvalue
	% No longitues here; thus half sum, with the sqrt(2) in there, is
        % OK; otherwise, full sum, with the sqrt(2) in there, would be
        % required, since we'd have cos(mphi) and sin(mphi) to deal with
	G=G+ones(size(V))*(E.^2)';
      end
      eval(sprintf('save %s F G',fnpl))
    end
  end
end

% Now plot them
clf
[ah,ha]=krijetem(subnum(2,2));
cols=[grey(7) ; 0 0 0 ; grey(5) ; grey(3)];
for ondex=1:length(TH)
  axes(ah(ondex))
  legsi{ondex}=sprintf('  %s = %i%s ','\Theta',TH(ondex),str2mat(176));

  for index=1:length(ZL)
    L=ZL(index);
    % Calculate N over A
    NA=(L+1)^2/4/pi;
   
    fnpl=sprintf('%s/SDSUMALL2-%i-%i-%i.mat',...
		 fullfile(getenv('IFILES'),'WIECZOREK'),...
		 TH(ondex),L,(L+1)^2);
    eval(sprintf('load %s',fnpl))
    
    % Verify that nothing exceeds N over A
    if max(F)>NA
      warning(sprintf('Weighted sum off N/A for L= %i and TH= %i',...
		      L,TH(ondex)))
    end
    % Verify that the unweighted G is indeed N over A
    if any(abs(G-NA)>1e-10)
      warning(sprintf(...
	  'Unweighted sum off N/A by %8.3f for L= %i and TH= %i',...
	  max(abs(G-NA)),L,TH(ondex)))
      disp('Use SDSUMK where this should have been fixed')
      disp('But verify what is wrong here later')
    end
    p(:,ondex)=plot(linspace(0,180,length(F)),F*100,...
		    'Color',cols(index,:),'LineW',1);
    hold on
    % At this point, 11-08-2005 I have tested the equivalency with SDSUMK
    % Comment this out or legend won't work
%    punw(ondex)=plot(linspace(0,180,length(G)),G*100,'LineW',1,...
%		     'Color',cols(index,:),'LineS','--'); 
  end
  xl(ondex)=xlabel(sprintf('colatitude %s','\theta'));
  yl(ondex)=ylabel('mean square error (%)');
  set(ah(ondex),'xlim',[0 180],'xtick',[0 TH(ondex) 180-TH(ondex) 180],...
		'ylim',[0 100],'ytick',[0:25:100],'xgrid','on','ygrid', ...
		'on')
  [bhh(ondex),thh(ondex)]=boxtex('ur',ah(ondex),legsi{ondex},12,[],0.85);
  deggies(ah(ondex),1)
  xlab=get(ah(ondex),'xtickl'); 
  xlab(end,:)=['    ']; xlab(1,:)=['    '];
  set(ah(ondex),'xtickl',xlab);
end

l=legend(legso,'Location','SouthEast');

% Cosmetics
longticks(ah)
nolabels(ha(3:4),2)
serre(ah(1:2),1/2,'across')
serre(ah(3:4),1/2,'across')
serre(ha(1:2),1/2,'down')
serre(ha(3:4),1/2,'down')

set([xl yl],'FontS',13)
set([ah],'FontS',12)

delete(xl(1:2))
delete(yl([2 4]))

fig2print(gcf,'portrait')
figdisp


