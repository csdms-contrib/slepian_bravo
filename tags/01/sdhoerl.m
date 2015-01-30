function sdhoerl(SN)
% SDHOERL(SN)
%
% Makes a mse/bias/variance plot that will enlighten us.
%
% Last modified by fjsimons-at-alum.mit.edu, 26.10.2005

defval('SN',4);
TH=10;
L=45;
NUM=100;

% Start the graphical layout
clf
[ah,ha]=krijetem(subnum(2,2));
fig2print(gcf,'portrait')
serre(ha(1:2),1/4,'down')
serre(ha(3:4),1/4,'down')
serre(ah(1:2),1/4,'across')
serre(ah(3:4),1/4,'across')

% Normalization factor for the signal strength
fax=(L+1)^2/4/pi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(ah(1))
% Calculate the contributions for damping levels
[eta,mseO,mseR,mseC,fO,fR,fC,varO,bs2O]=sdmse2(SN,NUM,TH,L);
difer(mseO-varO-bs2O)

eO=mseO/fax*100;
vO=varO/fax*100;
bO=bs2O/fax*100;

% Plot the averaged mse
h1(1)=plot(eta,eO,'k','LineW',1); hold on
% Plot the averaged variance
h1(2)=plot(eta,vO,'Color',grey,'LineW',1); 
h3(1)=plot([fO fO],[vO(find(eta==fO)) min(eO)],'k');

% Left Y-tickmarks and labels
% Widened range for the axis limits
yl1=1/SN*100+[-1 2.5];
% Where to put the ticks
ylt1=[yl1(1) vO(find(eta==fO)) min(eO) ...
      min(eO)+[yl1(2)-min(eO)]/2 yl1(2)];
% And round off the labels
yltm1=round(ylt1*10)/10;

set(ah(1),'ylim',yl1,'ycolor','k','box','off',...
	  'ytick',ylt1,'ytickl',yltm1,'xgrid','off',...
	  'ygrid','on','xtick',[0 fO 1],...
	  'xtickl',round([0 fO 1]*100)/100,...
	  'xgrid','on','ygrid','on')
xl(1)=xlabel(sprintf('damping parameter %s','\eta'));
yl(1)=ylabel(sprintf('%s-average mse, var, bias (%s)',...
			 '\Omega','%')); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(ah(2))
% Calculate the contributions for damping levels
[kay,mskO,mskR,mskC,kO,kR,kC,K,varkO,bsk2O]=sdmsk2(SN,[],TH,L);
difer(mskO-varkO-bsk2O)

% What is the range of kay's we really want to plot?
kr=50;
% Remember kay(1)=0
krg=[max(K-kr,1):min(K+kr,(L+1)^2)]+1;

tO=mskO/fax*100;
vkO=varkO/fax*100;
bkO=bsk2O/fax*100;

% Plot the averaged mse
h2(1)=plot(kay,tO); hold on
% Plot the averaged variance
h2(2)=plot(kay,vkO); 
h3(2)=plot([kO kO],[vkO(find(kay==kO)) min(tO)],'k');

% Left Y-tickmarks and labels
% Widened range for the axis limits
yl2=1/SN*100+[-1 2.5];
% Where to put the ticks
ylt2=[yl2(1) vkO(find(kay==kO)) min(tO) ...
      min(tO)+[yl2(2)-min(tO)]/2 yl2(2)];
% And round off the labels
yltm2=round(ylt2*10)/10;

set(ah(2),'ylim',yl2,'ycolor','k','box','off',...
	  'ytick',ylt2,'ytickl',yltm2,'xgrid','off',...
	  'ygrid','on','xlim',[kay(krg(1)) kay(krg(end))],...
	  'xtick',[kay(krg(1)) kO kay(krg(end))],...
	  'xdir','rev','xgrid','on','xgrid','on')
xl(2)=xlabel(sprintf('truncation rank k'));
yl(2)=ylabel(sprintf('%s-average mse, var, bias (%s)',...
			 '\Omega','%'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(ah(3))
% Now calculate curve for single theta
[mseTH,theta,varTH,bs2TH]=sderr(SN,TH,L,fO);
difer(mseTH-varTH-bs2TH)

half=[theta<=pi/2];
% Don't do square root no more
eTH=mseTH/fax*100;
varTH=varTH/fax*100;
bs2TH=bs2TH/fax*100;

h4(2)=plot(theta(half)*180/pi,varTH(half),'Color',grey,'LineW',1); 
hold on
h4(3)=plot(theta(half)*180/pi,eTH(half),'k','LineW',1); 
h4(1)=plot(theta(half)*180/pi,bs2TH(half),'k','LineW',0.5); 

xl(3)=xlabel(sprintf('colatitude %s','\theta'));
yl(3)=ylabel(sprintf('optimal mse, var, bias (%s)','%'));
set(ah(3),'xtick',[0 TH 45 90])
xlim([0 90])
set(ah(3),'ytick',unique(round([1/SN*100 25 50 75 100]*10)/10))
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(ah(4))
% Now calculate curve for single theta
keyboard
kO=320;
[mskTH,theta,varkTH,bsk2TH]=sderk(SN,TH,L,kO);
difer(mskTH-varkTH-bsk2TH)

half=[theta<=pi/2];
% Don't do square root no more
tkTH=mskTH/fax*100;
varkTH=varkTH/fax*100;
bsk2TH=bsk2TH/fax*100;

h5(2)=plot(theta(half)*180/pi,varkTH(half),'Color',grey,'LineW',1); 
hold on
h5(3)=plot(theta(half)*180/pi,tkTH(half),'k','LineW',1); 
h5(1)=plot(theta(half)*180/pi,bsk2TH(half),'k','LineW',0.5); 

% Cosmetics
set(h2,'marker','o','markers',1,'lines','none')
set(h2(1),'markerf','k','markere','k')
set(h2(2),'markerf',grey,'markere',grey)

set(ah,'box','on')

xl(4)=xlabel(sprintf('colatitude %s','\theta'));
yl(4)=ylabel(sprintf('optimal mse, var, bias (%s)','%'));
set(ah(4),'xtick',[0 TH 45 90])
xlim([0 90])
set(ah(4),'ytick',unique(round([1/SN*100 25 50 75 100]*10)/10))
grid on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set([xl yl],'FontS',13)
set(ah,'FontS',12)
deggies(ah(3:4),1)
longticks(ah)
[bh,th]=label(ah,'lr',12);
l(1)=legend(ah(1),'mse','var','bias^2','Location','NorthEast');
delete(yl([2 4]))

figdisp






