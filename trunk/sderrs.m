function sderrs(SN)
% sderrs(SN)
%
% Shows variance curves for a variety of SN ratios and dampings.
% Uses SDMSE2 to find the optimum eta and SDERR to calculate the MSE.
%
% Last modified by fjsimons-at-alum.mit.edu, 25.10.2005

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
  mineTH=10^9;
  % Figure out the minimum theta
  [eta,mseO,mseR,mseC,fO,fR,fC]=sdmse2(SN(ondex),[],TH,L);
  % Which damping levels do we plot?
  eta=[0 fO 1];
  legsi{ondex}=sprintf('N/S = %5.1f%s',1./SN(ondex)*100,'%');
  axes(ah(ondex))
  for index=1:length(eta)
    legso{index}=sprintf('%s = %5.2f','\eta',eta(index));
    [mseTH,theta]=sderr(SN(ondex),TH,L,eta(index),ntheta,NaN);

    half=[theta<=pi/2];
    % Don't do square root no more
    eTH=mseTH/fax*100;
    mineTH=min(mineTH,min(eTH));
    
    p(index,ondex)=plot(theta(half)*180/pi,eTH(half),...
		  'lines',lins{index},'color',cols{index});
    hold on
  end
  xl(ondex)=xlabel(sprintf('colatitude %s','\theta'));
  yl(ondex)=ylabel(sprintf('mse (%s)','%'));
  set(ah(ondex),'xtick',[0 TH 45 90])
  xlim([0 90])
  % ylim([mineTH 100]+[-1 1]*[100-mineTH]/15)
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

if length(SN)>1
  delete(xl(1:2))
  delete(yl([2 4]))
  nolabels(ah(1:2),1)
  deggies(ah(3:4),1)
end

set(p(:),'LineW',1)

movev(loh,.13)
longticks(ah(:))

figdisp











