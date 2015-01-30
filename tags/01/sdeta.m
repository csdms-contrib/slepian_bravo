function sdeta(SN)
% sdeta(SN)
%
% Plots damping vs variance curves from SDMSE.
%
% Last modified by fjsimons-at-alum.mit.edu, 21.10.2005

defval('SN',[4 6 8 10]);

TH=10;
L=45;
NUM=100;

% Start the graphical layout
clf
if length(SN)>1
  [ah,ha]=krijetem(subnum(length(SN)/2,2));
  fig2print(gcf,'portrait')
  serre(ha(1:2),1/4,'down')
  serre(ha(3:4),1/4,'down')
  % Special move for doubly annotated axes
  moveh(ha(1:2),-.01)
else
  ah=gca;
  fig2print(gcf,'portrait')
end

for index=1:length(SN)
  % Calculate the average variances and their predicted optima
  [eta,mseO,mseR,mseC,fO,fR,fC]=sdmse(SN(index),NUM,TH,L);

  legsi{index}=sprintf('N/S = %5.2f',1./SN(index));

  % Remember we are comparing bandlimited signal and error
  fax=(L+1)^2/4/pi;
  eO=sqrt(mseO/fax)*100;
  eR=sqrt(mseR/fax)*100;
  eC=sqrt(mseC/fax)*100;

  axes(ah(index))
  % Plot on two y-axes
  [ax(index,:),h1(index),h2(index)]=plotyy(eta,eO,eta,eR);

  % Left Y-tickmarks and labels
  % Widened range for the axis limits
  yl1=minmax(eO(3:end))+[-1 1]*range(eO(3:end))/15;
  % Where to put the ticks
  ylt1=[min(eO) min(eO)+[yl1(2)-min(eO)]/2 yl1(2)];
  % And round off the labels
  yltm1=round(ylt1*10)/10;

  % Right Y-tickmarks and labels
  yl2=minmax(eR(3:end))+[-1 1]*range(eR(3:end))/15;
  % Where to put the ticks
  ylt2=[min(eR) min(eR)+[yl2(2)-min(eR)]/2 yl2(2)];
  % And round off the labels
  yltm2=round(ylt2*100)/100;

  % Now set these tick marks and labels for the left axis
  set(ax(index,1),'ylim',yl1,'ycolor','k','box','off',...
		  'ytick',ylt1,'ytickl',yltm1,'xgrid','off',...
		  'ygrid','on','xtick',[0 fO 1],...
		  'xtickl',round([0 fO 1]*100)/100)
  xl(index)=xlabel(sprintf('damping parameter %s','\eta'));
  yl(index)=ylabel(sprintf('%s-average error-to-signal ratio (%s)',...
			   '\Omega','%')); 
  
  % And put the third axis on
  xx(index)=xtraxis(ah(index),fR,...
		    sprintf('%5.2f',fR));
  longticks(xx(index))
  
  % Now set these tick marks and labels for the right axis
  set(ax(index,2),'ylim',yl2,'xtick',[],'ycolor','k','box','off',...
		  'xaxisloc','top','ytick',ylt2,'ytickl',yltm2)
  
  yll(index)=ylabel(sprintf('%s-average error-to-signal ratio (%s)',...
			    'R','%'));
  hold on
  oO=plot([fO fO],minmax([yl1 yl2]),'-','LineW',0.5,'Color','k');
  oR=plot([fR fR],minmax([yl1 yl2]),'-','LineW',0.5,'Color',grey);

  [bh(index),th(index)]=boxtex('ur',ax(index,2),legsi{index},...
			       12,[],1.05,1.1);
  
end

set([h1 h2],'linew',1)
set(h1,'color','k')
set(h2,'color',grey)

set(xl,'FontS',13)
set([ah xx],'FontS',12)

if length(SN)>1
  delete(xl(1:2))
  delete(yl([2 4]))
  delete(yll([1 3]))
  moveh(yll([2 4]),.21)
  set(yll([2 4]),'rotation',270)
%  nolabels(ah(1:2),1)
end

longticks([ah(:) ;  ax(:)])

figdisp

