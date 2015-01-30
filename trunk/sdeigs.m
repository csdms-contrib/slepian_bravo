function sdeigs
% SDEIGS
%
% Simons and Dahlen (2006), Figure 8
% Plots eigenvalues of the eigenfunctions concentrated to the belt
% With the lightning bolt
% 
% Last modified by fjsimons-at-alum.mit.edu, 04/13/2007

warning('check out right tickmarks')

TH=[10 20 30 40]/2;
L=18;

clf
[ah,ha]=krijetem(subnum(2,2));

yls=[-0.1 1.1];

symbs={'o','x','s','+','v','*','^','d','<','>','p','h',...
       'o','s','v','d','p','h','^'};

more off
for index=1:length(TH)
  legsi{index}=sprintf('%s = %i%s','  \Theta',TH(index),str2mat(176));
  SN(index)=floor(L/180*TH(index));

  theta=0; phi=0; sord=2; srt='belt';
  
  [G,V,EM,GK,VK,NA,N(index)]=galpha(TH(index),L,sord,theta,phi,srt);

  % Total number to plot
  nuto=60;
  % Number of first ones to plot
  fones=11;
  % Minimum to plot in the second panel, after the break
  xmin=(L+1)^2-(nuto-fones);
  % Maximum possible value
  xmax=(L+1)^2;
  % What is the offset, effectively
  offs=xmin-fones-1;
 
  axes(ah(index))
  % Make sure these do not show through the other symbols
  plot(round([N(index) N(index)]-offs),yls,'k:')
  hold on
  plot([0 xmax],[0.5 0.5],'k:')
  plot([0 xmax],[0 0],'k:')
  plot([0 xmax],[1 1],'k:')

  % Plot the first couple of ones
  for ondi=1:fones
    theorder=abs(EM(ondi));
    if theorder>11
      colo='white'; else colo=grey;
    end
    p(ondi,index)=plot(ondi,V(ondi),symbs{theorder+1},'MarkerF',colo);
  end
  % Plot the ones from the second panel, after the break
  for newdex=xmin:xmax
    ondi=ondi+1;
    theorder=abs(EM(newdex));
    if theorder>11
      colo='white'; else colo=grey;
    end
    p(ondi,index)=plot(ondi,V(newdex),symbs{theorder+1},'MarkerF',colo);
  end
  
  % First panel tickmarks
  fpti=[1 10:10:100]; fpti=fpti(fpti<fones);
  
  % Second panel tickmarks
  lpti=sort([360:-10:xmin]);
  
  % Don't have the last one fall on the grid itself
  set(ah(index),'xlim',[0 xmax-offs+1],'ylim',yls,'ytick',[0:0.25:1],...
		'xgrid','off','ygrid','off')
  
  % Put the box on
  [bh(index),th(index)]=boxtex('ll',ah(index),legsi{index},12);
  
  % Put the break box on
  d(index)=fillbox([fones-1 fones yls],'w');

  % Set the new and adjusted tickmarks
  set(ah(index),'xtickl',[fpti lpti],'xtick',[fpti lpti-offs])
  drawnow
end

% Plot the order legend here
axes(ah(1))
xofs=22;
% A big box
fb=fillbox([xofs+2 xofs+36 0.75 -0.05],'w');
% The zeroth order
for ondi=1
  ypo=0+0.075*(ondi-1);
  pl(ondi,1)=plot(xofs+4,ypo,symbs{ondi},'MarkerF',grey);
  hold on
  tl(ondi,1)=text(xofs+7,ypo,sprintf('m =   %i',ondi-1),'FontS',8);
end
% Orders one through nine
for ondi=2:10
  ypo=0+0.075*(ondi-1);
  pl(ondi,1)=plot(xofs+4,ypo,symbs{ondi},'MarkerF',grey);
  hold on
  tl(ondi,1)=text(xofs+7,ypo,sprintf('m = %s %i','\pm',ondi-1),'FontS',8);
end
% Orders ten through eighteen
for ondi=11:19
  ypo=0+0.075*(ondi-1-9);
  pl(ondi,1)=plot(xofs+4+18,ypo,symbs{ondi},'MarkerF','white');
  hold on
  tl(ondi,1)=text(xofs+7+18,ypo,sprintf('m = %s %i','\pm',ondi-1),'FontS',8);
end

% Now make the plot beautiful
longticks(ah)
set([p(~~p(:)) ; pl(~~pl(:))],'MarkerS',4,'MarkerE','k')
axes(ha(1))
al(1)=ylabel('eigenvalue \lambda');
axes(ha(2))
al(2)=ylabel('eigenvalue \lambda');
axes(ha(2))
xl(1)=xlabel('rank \alpha');
axes(ha(4))
xl(2)=xlabel('rank \alpha');

nolabels(ha(3:4),2)
nolabels(ah(1:2),1)

serre(ah(1:2),1/2,'across')
serre(ah(3:4),1/2,'across')
serre(ha(1:2),1/2,'down')
serre(ha(3:4),1/2,'down')

for ind=1:4
  xx(ind)=xtraxis(ah(ind),round(N(ind)-offs),...
		  {sprintf('K = %i',round(N(ind)))});
% Not it  set(xx(ind),'box','on')
end
longticks(xx)

set([xl al],'FontS',13)
set([ah xx],'FontS',12)

fig2print(gcf,'portrait')
figdisp


