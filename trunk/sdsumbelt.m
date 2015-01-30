function sdsumbelt
% SDSUMBELT
%
% Simons and Dahlen (2006) Figure 6
% Sums all eigenfunctions on the LATITUDINAL BELT.
%
% Last modified by fjsimons-at-alum.mit.edu, 04/13/2007

TH=[10 20 30 40]/2;
L=18;

theta=linspace(0,pi,720);

% Plot and calculate at the same time
clf
[ah,ha]=krijetem(subnum(2,2));
cols=[grey(3) ; grey(6) ; 0 0 0 ];

% Overall legend
post={'1\rightarrow K' '1\rightarrow (L+1)^2' 'unweightedi'};

% Calculate and plot at the same time     
for index=1:length(TH)
  % Individual legend
  legsi{index}=sprintf('%s = %i%s    ','\Theta',TH(index),str2mat(176));
  
  % Calculate
  [G,V,EM,GK,VK,NA,N]=galpha(TH(index),L,2,theta,NaN,'belt');
  
  % Plot (un)weighted (un)truncated sums
  axes(ah(index))
  p(index,1)=plot(theta*180/pi,diag(G'*diag(V)*G),'Color',cols(1,:)); 
  hold on
  p(index,2)=plot(theta*180/pi,diag(GK'*diag(VK)*GK),'Color',cols(2,:));
  p(index,3)=plot(theta*180/pi,diag(G'*G),'Color',cols(3,:),'LineS','--');
  
  % A bit of cosmetics
  set(ah(index),'xtick',[0 TH(index) 90 180-TH(index) 180],'xgrid','on',...
		'ytick',[0 NA],'ygrid','off')
  xls=[0 90]; yls=[-1 35];
  set(ah(index),'xlim',xls,'ylim',yls)
  
  % Individual legend
  [bhh(index),thh(index)]=boxtex('ur',ah(index),legsi{index},13,0.85);
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
deggies(ah(1:4),1)

% This needs to be after deggies
set(ha(1:2),'ytickl',{'' ''})
set(ha(3:4),'ytickl',{'0' 'K/A'},'yaxisl','right')
set([yl xl],'FontS',15)
set(ah,'FontS',14)

set(p,'LineW',1)

serre(ah(1:2),1/2,'across')
serre(ah(3:4),1/2,'across')
serre(ha(1:2),1/2,'down')
serre(ha(3:4),1/2,'down')

% Now amend tick marks
% for index=1:4
%   xlab=get(ah(index),'xtickl'); 
%   xlab(1,:)=['    '];
%   xlab(end,:)=['    '];
%   set(ah(index),'xtickl',xlab);
% end
for index=1
  xlab=get(ah(index),'xtickl'); 
  xlab(1,:)=['   '];
  set(ah(index),'xtickl',xlab);
end

axes(ah(1))
loh(1)=legend(post,'Location','SouthEast');  
set(getkids(loh(1),2),'Color',cols(3,:),'LineS','--')
set(getkids(loh(1),5),'Color',cols(1,:))
set(getkids(loh(1),8),'Color',cols(2,:))

set(findobj('string','unweightedi'),'string','unweighted')

moveh(thh(:),5)

fig2print(gcf,'portrait')
figdisp

