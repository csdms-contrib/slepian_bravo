function sdkay2(SN)
% sdkay2(SN)
%
% Simons and Dahlen (2006), Figure 12
% Plots truncation vs variance curves from SDMSK2
% For mean square error, not its square root, and integration by hand.
%
% Last modified by fjsimons-at-alum.mit.edu, 09/23/2023

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

% Excessive verification
xver=0;
for index=1:length(SN)
  % Calculate the average variances and their predicted optima
  [kay,mskO,mskR,mskC,kO,kR,kC,K]=sdmsk2(SN(index),[],TH,L);
  if xver==1
    disp('Testing the difference between SDMSK and SDMSK2')
    [kay,mskO2,mskR2,mskC2]=sdmsk(SN(index),[],TH,L);
    difer(mskO-mskO2)
    difer(mskR-mskR2)
    difer(mskC-mskC2)
  end
  % Compare with the eta-damping approach
  [eta,mseO,mseR,mseC,fO,fR,fC]=sdmse2(SN(index),100,TH,L);
  if xver==1
    disp('Testing the difference between SDMSE and SDMSE2')
    [eta,mseO2,mseR2,mseC2]=sdmse(SN(index),100,TH,L);
    difer(mseO-mseO2)
    difer(mseR-mseR2)
    difer(mseR-mseR2)
  end
  
  % Check the consistency of undamped with untruncated
  difer(mseO(1)-mskO(end),[],0)
  difer(mseR(1)-mskR(end),[],0)
  
  % Write progressive legend
  legsi{index}=sprintf('N/S = %5.1f%s',1./SN(index)*100,'%');

  % Remember we are comparing bandlimited signal and error
  fax=(L+1)^2/4/pi;
  % Don't do square root no more
  tO=mskO/fax*100; eO=mseO/fax*100;
  tR=mskR/fax*100; eR=mseR/fax*100;
  tC=mskC/fax*100; eC=mseC/fax*100;

  % What is the range of kay's we really want to plot?
  kr=50;
  % Remember kay(1)=0
  krg=[max(K-kr,1):min(K+kr,(L+1)^2)]+1;
    
  axes(ah(index))
  % Plot on two y-axes
  [ax(index,:),h1(index),h2(index)]=...
      plotyy(kay(krg),tO(krg),kay(krg),tR(krg));
  % Predict the no truncation grey R-average
  difer(tR(end)-100/cos(TH*pi/180)/SN(index))
  
  % Left Y-tickmarks and labels
  % Widened range for the axis limits
  % Knowing that min(eO)<min(tO) and max(tO)>max(eO)
  yl1=1/SN(index)*100+[0 1]*5;
  % Where to put the tick marks, one at minimum, one half way
  ylt1=[yl1(1) min(tO) min(eO)+[yl1(2)-min(eO)]/2 yl1(2)];
  % And round off the labels
  yltm1=round(ylt1*10)/10;

  % Right Y-tickmarks and labels
  % Knowing that min(eR)<min(tR) and max(tR)>max(eR)
  % That's all good but want the minimum of grey to coincide with the
  % minum of the black one
  yl2=[1/SN(index)*100 ...
       min(tR)+(yl1(2)-min(tO))/(min(tO)-yl1(1))*(min(tR)-1/SN(index)*100)];
  % Where to put the tick marks, one at minimum, one half way
  ylt2=[yl2(1) min(tR) min(eR)+[yl2(2)-min(eR)]/2 yl2(2)];
  % And round off the labels
  yltm2=round(ylt2*10)/10;

  % Now set these tick marks and labels for the left axis
  set(ax(index,1),'ylim',yl1,'ycolor','k','box','off',...
		  'ytick',ylt1,'ytickl',yltm1,'xgrid','on',...
		  'ygrid','on',...
		  'xlim',[kay(krg(1)) kay(krg(end))],...
		  'xtick',[kay(krg(1)) K kay(krg(end))])

  xl(index)=xlabel(sprintf('truncation rank k'));
  yl(index)=ylabel(sprintf('%s-average mse (%s)','\Omega','%'));

  % And put the third axis on
  xx(index)=xtraxis(ah(index),kO,sprintf('%i',kO)); 
  set([ax(:) ; xx(:)],'xdir','rev')
  longticks(xx(index))

  % Now set these tick marks and labels for the right axis
  set(ax(index,2),'ylim',yl2,'xtick',[],'ycolor','k','box','off',...
		  'xaxisloc','top','ytick',ylt2,'ytickl',yltm2,...
		  'xlim',[kay(krg(1)) kay(krg(end))],...
		  'xtick',[])
  
  yll(index)=ylabel(sprintf('%s-average mse (%s)',...
			    'R','%'));
  hold on
  oR=plot([kR kR],minmax([yl1 yl2]),'-','LineW',1.5,'Color',grey);
  oO=plot([kO kO],minmax([yl1 yl2]),'-','LineW',0.5,'Color','k');
  difer(kO-kR)
  
  [bh(index),th(index)]=boxtex('ll',ax(index,2),legsi{index},...
			       12,[],1.05,1.1);
end

set([h1 h2],'marker','o','markers',2,'lines','none')
set(h1,'markerf','k','markere','k')
set(h2,'markerf',grey,'markere',grey)

set([xl(~~xl) yl(~~yl) yll(~~yll)],'FontS',13)
set([ah(:); ax(:) ; xx(:)],'FontS',12)

if length(SN)>1
  delete(xl(1:2))
  delete(yl([2 4]))
  delete(yll([1 3]))
  moveh(yll([2 4]),-.21*range(krg))
  set(yll([2 4]),'rotation',270)
end

longticks([ah(:) ;  ax(:)])

figdisp



