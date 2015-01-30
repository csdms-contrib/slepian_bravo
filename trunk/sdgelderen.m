function sdgelderen(TH,ZL)
% SDGELDEREN(TH,ZL)
%
% INPUT:
%
% TH      Colatitudinal radii of the concentration regions, vector
% ZL      Bandwidth, scalar or vector
%
% Makes a plot of the mean-square error in the spatial domain of a global
% function  recovered in the spherical harmonic basis from incomplete
% sampling due to a polar gap of axisymmetric antipodal polar caps. 
%
% See also SDEPSTH, SDSUMALL2, SDGELDEREN2
%
% Last modified by fjsimons-at-alum.mit.edu, 05/21/2009

defval('TH',[3 5 7 10])
defval('ZL',[45 90 180])

for index=1:length(TH)
  for ondex=1:length(ZL)
    % The bandwidth is adapted, and thus is nth in sdepsth
    L=ZL(ondex);
    legso{ondex}=sprintf('L = %i',L);
    % Get mean-square error, in percentage
    [e2th{index}{ondex},theta{ondex}]=sdepsth(TH(index),L);
    e2th{index}{ondex}=e2th{index}{ondex}*100;
  end
end

clf
[ah,ha]=krijetem(subnum(2,2));

cols=[grey(7) ; 0 0 0 ; grey(5) ; grey(3)];

for index=1:length(TH)
  axes(ah(index))
  legsi{index}=sprintf('  %s = %i%s','\Theta',TH(index),str2mat(176));
  for ondex=1:length(ZL)
    plot(theta{ondex}*180/pi,e2th{index}{ondex},'LineW',1,...
	 'Color',cols(ondex,:))	
    hold on
  end
  hold off
  xl(index)=xlabel(sprintf('colatitude %s','\theta'));
  yl(index)=ylabel('mean square error (%)');
  set(ah(index),'xlim',[0 180],'xtick',[0 TH(index) 180-TH(index) 180],...
		'ylim',[0 100],'ytick',[0:25:100],'xgrid','on','ygrid','on')
  [bh(index),th(index)]=boxtex('ur',ah(index),legsi{index},12,[],0.85);
  % Remove the second and the one but last labels
  deggies(ah(index),1)
  xlab=get(ah(index),'xtickl'); 
  xlab(end,:)=['    ']; xlab(1,:)=['    '];
  set(ah(index),'xtickl',xlab);
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

if length(xl)>2
  delete(xl(1:2))
  delete(yl([2 4]))
end

fig2print(gcf,'portrait')
figdisp
