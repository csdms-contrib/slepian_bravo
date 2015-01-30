function sdgrunval
% SDGRUNVAL
%
% Simons & Dahlen (2006), Figure 9
% Eigenfunctions of the Grunbaum operator for the double polar cap
%
% Last modified by fjsimons-at-alum.mit.edu, 05.05.2006

TH=[10 20 30 40]/2;
L=18;

clf
[ah,ha]=krijetem(subnum(2,2));

% Calculations
for index=1:length(TH)
  axes(ah(index))
  legsi{index}=sprintf(' %s = %i%s ','  \Theta',TH(index),str2mat(176));
  for ondex=0:L
    [E,V,th,C,T]=grunbaum2(TH(index),L,ondex,0);
    p(ondex+1,index)=plot((ondex+1):(L+1),V-50*ondex,'Marker',symbol(1,1));
    pp(ondex+1)=V(end)-50*ondex;
    hold on
    drawnow
  end
  ylim([-1300 50])
  [bh(index),th(index)]=boxtex('ul',ah(index),legsi{index},12);
end

% Cosmetics
longticks(ah)

axes(ha(1))
al(1)=ylabel('eigenvalue \chi');
axes(ha(2))
al(2)=ylabel('eigenvalue \chi');
axes(ha(2))
xl(1)=xlabel('rank + order \alpha + m');
axes(ha(4))
xl(2)=xlabel('rank + order \alpha + m');

nolabels(ha(3:4),2)
nolabels(ah(1:2),1)

serre(ah(1:2),1/2,'across')
serre(ah(3:4),1/2,'across')
serre(ha(1:2),1/2,'down')
serre(ha(3:4),1/2,'down')

set([xl al],'FontS',15)
set([ah],'FontS',14)

set(ah,'xgrid','on','ygrid','on','xtick',[1 5 10 15 20],'xlim',[0 22])
set(p,'color','k')
set(p,'MarkerS',4,'MarkerF',grey,'MarkerE','k')

axes(ah(4))
for index=0:2:L
  ti(index+1)=text(20,pp(index+1),sprintf('%i',index));
end

fig2print(gcf,'portrait')
figdisp

