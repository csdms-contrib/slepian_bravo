function sdupcont(wirt)
% SDUPCONT(wirt)
%
% Shows the effect of upward/downward continuation of the Slepian
% functions 
%
% Last modified by fjsimons-at-alum.mit.edu, 04/13/2007

TH=10;
L=18;
sord=2;
yls=[-0.09 1.1];

clf
[ah,ha]=krijetem(subnum(2,2));
fig2print(gcf,'portrait')
serre(ha(1:2),1/4,'down')
serre(ha(3:4),1/4,'down')
serre(ah(1:2),1/2,'across')
serre(ah(3:4),1/2,'across')

a=[-1 -2/3 -1/3 0 1/3 2/3 1];

theta=linspace(0,180,720);

fnpl=sprintf('%s/SDUPCONT-%i-%i-%i.mat',...
	     fullfile(getenv('IFILES'),'SDWCAP'),TH,L,sord);
if exist(fnpl,'file')==2 & 1==3
  eval(sprintf('load %s',fnpl))
  disp(sprintf('%s loaded by SDUPCONT',fnpl))
else
  for index=1:length(a)
    % Do not sort, want to keep original sorting for comparison with a=0
    [G{index},V{index},EM{index},jk2,jk3,NA,N]=...
	galpha(TH,L,2,[],NaN,'local',a(index),1);
  end
  eval(sprintf('save %s G V EM',fnpl))
end

% Which of the eigenfunctions are we plotting?
% wirt=max(1,round(rand(1)*(L+1)^2));
defval('wirt',1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xmax=(L+1)^2;
% Eigenvalue behavior
axes(ah(1))
for index=1:4
  % Make a global sort on the BELT now
  [vv,jj]=sort(1-V{index},'descend');
  e(index)=plot(vv,'-o');
  K(index)=sum(1-V{index});
  hold on
  xlim([310 xmax])
%  wirt=jj(351);
%  1-V{index}(wirt)
end
plot([0 xmax],[0.5 0.5],'k:')
plot([0 xmax],[0 0],'k:')
plot([0 xmax],[1 1],'k:')
plot(round([K(4) K(4)]),yls,'k:')
% Cannot plot wirt as the eigenvalues plotted are now globally sorted 
hold off

axes(ah(2))
for index=4:7
  % Make a global sort on the BELT now
  [vv,jj]=sort(1-V{index},'descend');
  f(index-3)=plot(vv,'-o');
  K(index)=sum(1-V{index});
  hold on
  xlim([310 xmax])
end
top(f(1),gca)
plot([0 xmax],[0.5 0.5],'k:')
plot([0 xmax],[0 0],'k:')
plot([0 xmax],[1 1],'k:')
plot(round([K(4) K(4)]),yls,'k:')
hold off

axes(ah(3))
for index=1:4
  p(index)=plot(theta,G{index}(wirt,:));
  hold on
end

axes(ah(4))
for index=4:7
  u(index-3)=plot(theta,G{index}(wirt,:));
  hold on
end

% Cosmetics
set(p(end),'linew',1)
set(p(1:end-1),'linew',0.5)
set([p(1) e(1)],'color',grey(7))
set([p(2) e(2)],'color',grey(5))
set([p(3) e(3)],'color',grey(2))
set([p(end) e(end)],'color','k')

set(u(1),'linew',1)
set(u(2:end),'linew',0.5)
set([u(4) f(4)],'color',grey(7))
set([u(3) f(3)],'color',grey(5))
set([u(2) f(2)],'color',grey(2))
set([u(1) f(1)],'color','k')

set(e,'MarkerS',2)
set(e(1),'Markere',grey(7),'MarkerF',grey(7))
set(e(2),'MarkerE',grey(5),'MarkerF',grey(5))
set(e(3),'MarkerE',grey(2),'MarkerF',grey(2))
set(e(end),'MarkerE','k','MarkerF','k')

set(f,'MarkerS',2)
set(f(4),'MarkerE',grey(7),'MarkerF',grey(7))
set(f(3),'MarkerE',grey(5),'MarkerF',grey(5))
set(f(2),'Markere',grey(2),'MarkerF',grey(2))
set(f(1),'MarkerE','k','MarkerF','k')

set(f(1),'MarkerE','k','MarkerF','k')
set(ah(1:2),'ylim',yls,'ytick',[0:0.25:1],'ygrid','off')
set(ah(3:4),'xlim',[0 180],'xtick',[TH 90 180-TH],'xgrid','on')
set(ah(3:4),'ylim',[-1.25 4],'ytick',[0 2 4],'ygrid','on')

deggies(ah(3:4),1)
nolabels(ha(3:4),2)
longticks(ah)
[bh,th]=label(ah,'ll',12);

axes(ah(1))
yl(1)=ylabel('eigenvalue \lambda');
xl(1)=xlabel('rank \alpha');
tl(1)=title('downward continuation');
l(1)=legend('a =  1','a = 2/3','a = 1/3','a =  0',...
	    'location','South');

axes(ah(2))
xl(2)=xlabel('rank \alpha');
tl(2)=title('upward continuation');

axes(ah(3))
yl(3)=ylabel(sprintf('worst m = %i eigenfunction',EM{1}(wirt)));
xl(3)=xlabel('colatitude \theta');
l(2)=legend('a =  1','a = 2/3','a = 1/3','a =  0',...
	    'location','North');

axes(ah(4))
xl(4)=xlabel('colatitude \theta');
figdisp

set([tl],'FontS',15)
set([xl(~~xl) yl(~~yl)],'FontS',13)
set([ah],'FontS',12)




 
