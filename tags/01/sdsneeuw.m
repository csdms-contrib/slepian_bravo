function sdsneeuw(TH,L,bas)
% SDSNEEUW(TH,L,bas)
%
% INPUT:
%
% TH      Colatitudinal radii of the concentration regions, vector
% L       Bandwidth, scalar
% bas     Basis, 'sh' (default) or 'slepian'
%
% Makes a plot of the mean-square error wedge of spherical harmonic
% coefficients of a global function recovered from incomplete sampling
% under the wrong minimization criterion, thus utter garbage, really.
%
% See also SDEPSLM, SPIE2009_2
%
% Last modified by fjsimons-at-alum.mit.edu, 05/21/2009

defval('L',90)
defval('TH',[3 5 7 10])
defval('bas','sh')

for index=1:length(TH)
  % Get mean-square spectral error, in percentage
  e2lm{index}=sdepslm(TH(index),L,bas)*100;
end

clf
[ah,ha]=krijetem(subnum(2,4));

% Don't plot all the orders, just the first mor
mor=20;

colh=flipud(gray(20));
cols='gray(20)';

for index=1:4
  axes(ah(index))
  legsi{index}=sprintf(' %s = %i%s','\Theta',TH(index),str2mat(176));
  imagefnan([0 0],[mor L],e2lm{index}(:,1:mor+1),colh, ...
	    minmax(e2lm{index}(:,1:mor+1)))
  axis ij
  xl(index)=xlabel('order |m|');
  switch bas
    case 'sh'
     yl(index)=ylabel('degree l');
   case 'slepian'
    yl(index)=ylabel('rank + order \alpha + |m|');
  end
  [bh(index)]=title(legsi{index});
  
  % Add a line demarcating the physically possible
  hold on
  plaus(index)=plot([0 mor+1],[-1 mor],'k-','LineW',0.5);
  % And the turning point - where does this come from again?
  pturn(index)=plot((sin(TH(index)*pi/180)*([-1 L+1]+1/2)),[-1 L+1],'k-',...
		    'LineW',0.5);
  hold off
end

% Cosmetics
longticks(ah)
nolabels(ah(2:4),2)
serre(ah(1:4),2.5,'across')

set([xl yl],'FontS',11)
set([ah],'FontS',10)

moveh(ah,-.15)
movev(ah,-.15)

delete(yl([2:4]))
delete(ah(5:8))

caxcon=[0 100];
caxoc=[0 100];
[cb1,xcb1]=addcb('vert',caxcon,caxoc,cols,20);

set(cb1,'YaxisL','r')
set(xcb1,'string','mean square error (%)')
moveh(cb1,-0.02)

delete(xl([1 2 4]))
moveh(xl(3),-12)

fig2print(gcf,'portrait')

if strcmp(bas,'sh')
  figdisp
elseif strcmp(bas,'slepian')
  figdisp('sdsneeuw2')
end
