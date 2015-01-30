function sdsumall
% SDSUMALL
%
% Simons & Dahlen and Wieczorek (2005)
%
% Sums all eigenfunctions on the BELT.
%
% First calculate a bunch of cumulative sums
%
% For large bandwidths and small angles, see SDSUMALL2.
%
% Last modified by fjsimons-at-alum.mit.edu, 09.08.2005
% Stupid program, does double work. Fix later. See SDSUMK.
% ALSO, WE JUST ADAPT IT LINE BY LINE TO OUR NEEDS. DUMB.

TH=[10 20 30 40]/2;
L=18;
nth=128;

verify=0;

disp('VERIFY!')

disp(' ')
disp('Output from SDSUMALL')
disp(' ')

for ondex=1:length(TH)
  fnpl=sprintf('%s/SDSUMCUM-%i-%i-%i.mat',...
	       fullfile(getenv('IFILES'),'WIECZOREK'),...
	       TH(ondex),L,(L+1)^2);
  % Calculate the Shannon number for the BELT
  N(ondex)=round((L+1)^2*(1-(1-cos(TH(ondex)*pi/180))/2*2));
  N(ondex)=round((L+1)^2*((1-cos(TH(ondex)*pi/180))/2*2))
  
  % This is a fix on what is fast becoming an ugly program.
  % Since N may split a pair of degenerate eigenvalues in the middle
  if TH(ondex)==10; N(ondex)=N(ondex)+1; end
  
  NA=(L+1)^2/4/pi;

  disp(sprintf('TH = %2.2i ; Shannon number is %3.3i ; N/A = %8.3f',...
	       TH(ondex),N(ondex),NA))

  if exist(fnpl,'file')~=2  | 1==1
    % Figure out order to sum them in; the 2 is for the BELT
%    [lrnk,mrnk,lval]=sdelm(TH(ondex),L,2,1);
    [lrnk,mrnk,lval]=sdelm(TH(ondex),L,1,1);

    % Must have all the pairs of orders to be longitudinally independent
    if mrnk(N(ondex))~=0 & ...  % Not zero
	  lrnk(N(ondex))==lrnk(N(ondex)+1) & ... % Same rank
	  [mrnk(N(ondex))==-mrnk(N(ondex)+1)] % Just switched signs
      N(ondex)=N(ondex)+1;
    end
    
    % Calculate cumulative sums weighted by eigenvalue!
    F=0;
    for index=1:length(lrnk)
      % Note this is not very efficient since all the nonzero orders are
      % simpy done twice!
      % This is at different longitudes, of which 8 are plotted
      % Its over all +/- orders,  the cos and sin must cancel 
      [E,Vg,th,C,T,V]=grunbaum2(TH(ondex),L,mrnk(index),nth,2);
      % Watch for this eigenvalue! It's for the BELT
      % And, the factor two is in there
%      F=F+(E{lrnk(index)}).^2.*(1-V(lrnk(index)));
      F=F+(E{lrnk(index)}).^2.*(V(lrnk(index)));
      if index==length(lrnk) | index==N(ondex)
	fnpl=sprintf('%s/SDSUMCUM-%i-%i-%i.mat',...
		     fullfile(getenv('IFILES'),'WIECZOREK'),...
		     TH(ondex),L,index);
	eval(sprintf('save %s F',fnpl))
      end
    end
  end
end

if verify==1
  disp(' ')
  disp('Output from SDSUMK')
  disp(' ')
end

% Now plot them
clf
[ah,ha]=krijetem(subnum(2,2));
xls=[0 180]; yls=[-1 35];
cols=[grey(3) ; grey(6) ; 0 0 0 ];

for ondex=1:length(TH)
  % WATCH OUT, LEGENDS ARE LATER MOVEH'D
  legsi{ondex}=sprintf(' %s = %i%s  ','\Theta',TH(ondex),str2mat(176));
  % Sum of how many terms are plotted? 
  pott=[(L+1)^2 N(ondex)];
  post{ondex}=...
     {'1\rightarrow K' '1\rightarrow (L+1)^2' 'unweightedi'};
  
  axes(ah(ondex))
  
  % To verify both complete and Shannon-number sums
  if verify==1
    [F2,G2,N2,NA2,th2]=sdsumk(TH(ondex),L,nth,2);
  end

  for index=1:length(pott)
    fnpl=sprintf('%s/SDSUMCUM-%i-%i-%i.mat',...
		 fullfile(getenv('IFILES'),'WIECZOREK'),...
		 TH(ondex),L,pott(index));
    eval(sprintf('load %s',fnpl))
    
    % We plot 8 longitudes, which, however, should all overlap
    longtol=sum(sum(abs(diff(F,1,2)),2))/prod(size(F));
    if longtol>100*eps
      disp(sprintf('Longitudinal error tolerance exceeded ; error is %8.3f',... 
		   longtol))
    end
    noteight=1;
    p{index}(:,ondex)=plot(linspace(0,180,size(F,1)),...
	 F(:,round(linspace(1,size(F,2),noteight))),...
	 'Color',cols(index,:));
    hold on

    if verify==1
      % Verify Shannon number sum
      if index==2
	plot(th2,F2,'b-')
	verif=sum(F(:,round(rand*size(F,2)))-F2);
	if verif>1e-10
	  warning(sprintf(...
	      'SDSUMALL and SDSUMK do not agree on K, %8.3f',verif))
	end
      elseif index==1
      % Verify complete number sum
	plot(th2,G2,'g-')
	verif=sum(F(:,round(rand*size(F,2)))-G2);
	if verif>1e-10
	  warning(sprintf(...
	      'SDSUMALL and SDSUMK do not agree on full sum, %8.3f',verif))
	end
      end
    end
  end  
  
 % First define general tickmarks for all panels
  set(ah(ondex),'xtick',[0 TH(ondex) 90 180-TH(ondex) 180],'xgrid','on',...
		'ytick',[0 NA],'ygrid','off')
  set(ah,'xlim',xls,'ylim',yls)

  G=0;
  for jndex=0:L
    [E2,Vg2,th,C2,T2,V2]=grunbaum2(TH(ondex),L,jndex,nth);
    % Remember here too, the factor of sqrt(2) was already in there
    G=G+sum(E2.^2,2);
  end
  punw(ondex)=plot(th,G,'LineW',1,'Color',cols(3,:),'LineS','--');
  % END ADD UNWEIGHTED STUFF
  [bhh(ondex),thh(ondex)]=boxtex('ur',ah(ondex),legsi{ondex},12,0.85);
end
longticks(ah)

axes(ha(1))
yl(1)=ylabel('cumulative energy');
axes(ha(2))
yl(2)=ylabel('cumulative energy');
axes(ah(3))
xl(1)=xlabel('colatitude \theta');
axes(ah(4))
xl(2)=xlabel('colatitude \theta');
deggies(ah(2:4),1)

deggies(ah(1),1)
% This needs to be after deggies
set(ha(1:2),'ytickl',{'' ''})
set(ha(3:4),'ytickl',{'0' 'N/A'},'yaxisl','right')
set([yl xl],'FontS',13)
set(ah,'FontS',12)

set([p{:}],'LineW',1)
serre(ah(1:2),1/2,'across')
serre(ah(3:4),1/2,'across')
serre(ha(1:2),1/2,'down')
serre(ha(3:4),1/2,'down')

for index=1
  axes(ah(index))
  loh(index)=legend(post{index},'Location','SouthEast');  
  set(getkids(loh(index),2),'Color',cols(3,:),'LineS','--')
  set(getkids(loh(index),5),'Color',cols(1,:))
  set(getkids(loh(index),8),'Color',cols(2,:))
end
%movev(loh,-.070)
moveh(thh(:),5)

% Now amend tick marks
for index=1:4
  xlab=get(ah(index),'xtickl'); 
  xlab(1,:)=['    '];
  xlab(end,:)=['    '];
  set(ah(index),'xtickl',xlab);
end

set(findobj('string','unweightedi'),'string','unweighted')

fig2print(gcf,'portrait')
figdisp


