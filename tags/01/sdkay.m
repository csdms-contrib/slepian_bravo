function sdkay(SN)
% sdkay(SN)
%
% Plots truncation vs variance curves from SDMSK (and SDMSE)
%
% See also SDKAY2
%
% Last modified by fjsimons-at-alum.mit.edu, 24.10.2005

defval('SN',[4 6 8 10]);

TH=10;
L=45;

% Start the graphical layout
clf
if length(SN)>1
  [ah,ha]=krijetem(subnum(length(SN)/2,2));
  fig2print(gcf,'portrait')
  serre(ha(1:2),1/4,'down')
  serre(ha(3:4),1/4,'down')
  moveh(ha(1:2),-.01)
else
  ah=gca;
  fig2print(gcf,'portrait')
end

for index=1:length(SN)
  % Calculate the average variances and their predicted optima
  % This integrated the functions for real:
  [kay,mskO,mskR,mskC,kO,kR,kC,K]=sdmsk(SN(index),[],TH,L);
  % Compare with the eta-damping approach
  % This integrated the functions for real:
  [eta,mseO,mseR,mseC,fO,fR,fC]=sdmse(SN(index),100,TH,L);

  % Check the consistency of undamped with untruncted
  difer(mseO(1)-mskO(end),[],0)
  difer(mseR(1)-mskR(end),[],0)
  
  legsi{index}=sprintf('N/S = %5.2f',1./SN(index));

  % Remember we are comparing bandlimited signal and error
  fax=(L+1)^2/4/pi;
  tO=sqrt(mskO/fax)*100; eO=sqrt(mseO/fax)*100;
  tR=sqrt(mskR/fax)*100; eR=sqrt(mseR/fax)*100;
  tC=sqrt(mskC/fax)*100; eC=sqrt(mseC/fax)*100;

  % What is the range of kay's we really want to plot?
  kr=50;
  % Remember kay(1)=0
  krg=[max(K-kr,1):min(K+kr,(L+1)^2)]+1;
    
  axes(ah(index))
  % Plot on two y-axes
  [ax(index,:),h1(index),h2(index)]=...
      plotyy(kay(krg),tO(krg),kay(krg),tR(krg));

  % Left Y-tickmarks and labels
  % Widened range for the axis limits
  % Knowing that min(eO)<min(tO) and max(tO)>max(eO)
  yl1=[min(eO) max(tO(krg(1:end-5)))+range(tO(krg(1:end-5)))/25];
  % Where to put the tick marks, one at minimum, one half way
  ylt1=[min(tO) min(eO)+[yl1(2)-min(eO)]/2 yl1(2)];
  % And round off the labels
  yltm1=round(ylt1*10)/10;

  % Right Y-tickmarks and labels
  % Widened range for the axis limits
  % Knowing that min(eR)<min(tR) and max(tR)>max(eR)
  % That's all good but want the minimum of grey to coincide with the
  % minum of the black one
  yl2=[min(eR) ...
       min(tR)+(yl1(2)-ylt1(1))/(ylt1(1)-yl1(1))*(min(tR)-min(eR))];

  % Where to put the tick marks, one at minimum, one half way
  ylt2=[min(tR) min(eR)+[yl2(2)-min(eR)]/2 yl2(2)];
  % And round off the labels
  yltm2=round(ylt2*100)/100;

  set(ax(index,1),'ylim',yl1,'ycolor','k','box','off',...
		  'ytick',ylt1,'ytickl',yltm1,'xgrid','on',...
		  'ygrid','on',...
		  'xlim',[kay(krg(1)) kay(krg(end))],...
		  'xtick',[kay(krg(1)) K kay(krg(end))])
  
  xl(index)=xlabel(sprintf('truncation level'));
  yl(index)=ylabel(sprintf('%s-average error-to-signal ratio (%s)',...
			   '\Omega','%'));
  xx(index)=xtraxis(ah(index),kO,...
		    sprintf('%i',kO));
  set([ax(:) ; xx(:)],'xdir','rev')
  longticks(xx(index))

  set(ax(index,2),'ylim',yl2,'xtick',[],'ycolor','k','box','off',...
		  'xaxisloc','top','ytick',ylt2,'ytickl',yltm2,...
		  'xlim',[kay(krg(1)) kay(krg(end))],...
		  'xtick',[])
  
  yll(index)=ylabel(sprintf('%s-average error-to-signal ratio (%s)',...
			    'R','%'));
  hold on
  oR=plot([kR kR],minmax([yl1 yl2]),'-','LineW',1.5,'Color',grey);
  oO=plot([kO kO],minmax([yl1 yl2]),'-','LineW',0.5,'Color','k');
  difer(kO-kR)
  
  [bh(index),th(index)]=boxtex('ul',ax(index,2),legsi{index},...
			       12,[],1.05,1.1);
end

set([h1 h2],'marker','o','markers',2,'lines','none')
set(h1,'markerf','k','markere','k')
set(h2,'markerf',grey,'markere','k')

set(xl,'FontS',13)
set([ah xx],'FontS',12)

if length(SN)>1
  delete(xl(1:2))
  delete(yl([2 4]))
  delete(yll([1 3]))
  moveh(yll([2 4]),-.21*range(krg))
  set(yll([2 4]),'rotation',270)
end

longticks([ah(:) ;  ax(:)])

figdisp





