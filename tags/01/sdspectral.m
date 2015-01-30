function sdspectral 
% SDSPECTRAL
%
% Simons & Dahlen (2005)
% Plots spectral functions for DOUBLE spherical polar cap
% for various angular orders with eigenvalue grey shading. 
%
% Last modified by fjsimons-at-alum.mit.edu, 22.07.2005

TH=30;
nth=128;
L=18; SN=6;
M=4;
Lnyq=nth-1;

clf
[ah,ha]=krijetem(subnum(M+1,SN));
fig2print(gcf,'landscape')

yls=[-120 10];

for m=0:M
  [E1,V1{m+1},th1,C1]=sdwcap2(TH,L,m,nth);
  [E2,V2{m+1},th2,C2]=sdwcapt2(TH,L,m,nth);

  E1m(m+1)=max(abs(E1(:)));
  E2m(m+1)=max(abs(E2(:)));
  narm2=1;

  for ondex=1:SN
    axes(ah(m*SN+ondex))
    warning off
    spec=decibel(narm2.*C2(:,ondex).^2);
    warning off
    spec(isinf(spec))=NaN;
    pg(m+1,ondex)=plot([m:Lnyq],spec,'o','Color',grey(6),...
		       'MarkerE',grey(6),'MarkerF',grey(6),'markers',1);
    hold on
    pm(m+1,ondex)=plot([m:L],spec(1:L-m+1),'o','Color','k',...
		       'MarkerE','k','MarkerF','k','markers',1);
    hold on
    plot([m m],yls,'k:')
    hold off
  end
end

nolabels(ah(1:end-SN),1)
nolabels(ha(M+2:end),2)
set(ah,'xlim',[0 Lnyq],'xtick',[0 L Lnyq],'xgrid','on',...
       'ygrid','on','ylim',yls)

% Eigenvalue labels
nf=9;
set(ah,'Fonts',nf)

for m=0:M
  for ondex=1:SN
    axes(ah(m*SN+ondex))
        t{m+1,ondex}=sprintf('%s = %9.6f','\lambda',V2{m+1}(ondex));
	if V2{m+1}(ondex)<0.975
	  lox='lr';
	else
	  lox='ur';
	end
	[bh(m+1,ondex),th(m+1,ondex)]=...
	    boxtex(lox,ah(m*SN+ondex),t{m+1,ondex},nf-1,[],0.85,1);
  end  
end


% Cosmetics
set(th,'FontS',nf-1)

longticks(ah,1/2)

for index=SN*M+1:SN*(M+1)
  axes(ah(index))
  xl1(index)=xlabel(sprintf('degree l'));
end

for index=1:SN
  axes(ah(index))
  tlb(index)=title(sprintf('%s = %i','\alpha',index));
end

set([pm(:) ; pg(:)],'LineW',1)

for index=1:M+1
  serre(ah([1:SN]+(index-1)*SN),1/3)
  axes(ha(index))
  ylb(index)=ylabel(sprintf('m = %i',index-1));
end

% This only for landscape
for index=1:SN
  serre(ha([1:M+1]+(index-1)*(M+1)),1,'down')
end

% In interactive mode this returns the right result, not on its own ?!?!?!

shrink(ah,1,1.2)

axes(ha(M+1))
tt=text(-60,-170,'dB','FontS',nf);

figdisp 


