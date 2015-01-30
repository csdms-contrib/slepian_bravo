function sderks(SN)
% sderrs(SN)
%
% Shows variance curves for a variety of SN ratios and truncations.
% Uses SDMSK2 to find the optimum kay and SDERK to calculate the MSE.
%
% Last modified by fjsimons-at-alum.mit.edu, 10.02.2006

defval('SN',[4 6 8 10]);

TH=10;
L=45;
ntheta=720;

fax=(L+1)^2/4/pi;

lins={'-','-','--'};
cols={grey,'k',grey};

clf
if length(SN)>1
  [ah,ha]=krijetem(subnum(length(SN)/2,2));
  fig2print(gcf,'portrait')
  serre(ah(1:2),1/3,'across')
  serre(ah(3:4),1/3,'across')
  serre(ha(1:2),1/2,'down')
  serre(ha(3:4),1/2,'down')
else
  ah=gca;
  fig2print(gcf,'portrait')
end

for ondex=1:length(SN)
  minkTH=10^9;
  % Figure out the minimum kay
  [kay,mskO,mskR,mskC,kO,kR,kC,K]=sdmsk2(SN(ondex),[],TH,L);
  % Which truncation levels do we plot?
  kay=[kay(end) kO K];
  legsi{ondex}=sprintf('N/S = %5.1f%s',1./SN(ondex)*100,'%');
  axes(ah(ondex))
  for index=1:length(kay)
    legso{index}=sprintf('%s = %i','J',kay(index));
    [mskTH,theta]=sderk(SN(ondex),TH,L,kay(index),ntheta,NaN);
    half=[theta<=pi/2];
    % Don't do square root no more
    kTH=mskTH/fax*100;
    minkTH=min(minkTH,min(kTH));
    
    p(index,ondex)=plot(theta(half)*180/pi,kTH(half),...
		  'lines',lins{index},'color',cols{index});
    hold on
  end
  xl(ondex)=xlabel(sprintf('colatitude %s','\theta'));
  yl(ondex)=ylabel(sprintf('mse (%s)','%'));
  set(ah(ondex),'xtick',[0 TH 45 90])
  xlim([0 90])
  yl1=[1/SN(ondex)*100 100];
  ylim(yl1+[-1 0]*range(yl1)/15)
  set(ah(ondex),'ytick',unique(round([yl1(1) 25 50 75 yl1(2)]*10)/10))

  [bh(ondex),th(ondex)]=boxtex('ur',ah(ondex),legsi{ondex},...
			       12,[],1.05,1.1);
  loh(ondex)=legend(ah(ondex),legso,'Location','SouthEast');  
  grid on
end
hold off

set([xl yl],'FontS',13)
set(ah,'FontS',12)
movev(loh,0.13)

if length(SN)>1
  delete(xl(1:2))
  delete(yl([2 4]))
  nolabels(ah(1:2),1)
  deggies(ah(3:4),1)
end

set(p(:),'LineW',1)

longticks(ah(:))

figdisp












