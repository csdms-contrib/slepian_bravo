function sdcmb4
% SDCMB4
%
% Simons & Dahlen (2005)
% Last four of the belt functions for a variety of orders
%
% Last modified by fjsimons-at-alum.mit.edu, 04/13/2007

TH=30;
L=18;
nth=32;
nlon=2*nth-1;

% The orders to display
m=[0 1 2];
% The number of functions for each order
numf=4;

% Collect positions for later overlap
[ah2,ha2]=krijetem(subnum(length(m),numf));
for index=1:length(m)*numf
  lox(index,:)=get(ah2(index),'position');
end
clf
[ah,ha]=krijetem(subnum(length(m),numf));

% Find all the eigenfunctions with all orders
theta=linspace(0,pi,nth);
phi=linspace(0,2*pi,nlon);
[G,V,EM,GK,VK,NA,N]=galpha(TH,L,2,theta,phi,'local');

% Loop over the orders and collect printed eigenvalues
numo2=0; Valo=[];
for ord=m
  % CAP FUNCTIONS
  for index=1:numf
    % Get the appropriate eigenfunction for the belt
    E=reshape(G(numo2+index,:),[nth nlon]);
    % Make sure the sign is right
    if E(1,round(end/2))<0
      E=-E;
    end
        % and some special flippings to coincide visually with the rest
    if (ord==1&index==1)|(ord==1&index==2)|(ord==2&index==3)|(ord==2&index==4)
      E=-E;
    end
    % Make a running index
    rindex=index+ord*numf;
    % Plot on the sphere
    axes(ah(rindex))
    plotonsphere(E,0.15*(1+[ord~=0]))
    view(90,20)
    % Set appropriate title string
    Valo=[Valo 1-V(numo2+index)];
    % Are we doing all right?
    if EM(numo2+index)~=ord
      error('You sure you got the numbering right?')
    end
    % Collect the range of values
    EMAX(rindex)=max(max(abs(E)));
  end
  % Cumulative number of functions per order to count forwards
  numo2=numo2+2*(L-abs(ord))+1;
end

% Movem all down
%for index=1:numf
%  serre(ha([1:length(m)]+(index-1)*length(m)),1/3,'down')
%end

% Put the titles on
for index=1:length(m)
  [lx(index),loc]=laxis(ha(index));
  ylb(index)=text(loc(1),loc(2),sprintf('m = %i',m(index)));
  set(ylb(index),'Rotation',90)
end
for index=1:numf
  [lxa(index),loc]=laxis(ah(index),0);
  % Set title string
  tits{index}=sprintf('L-m%+i',2-index);
  tt(index)=title(sprintf('%s = %s','\alpha',tits{index}));
end

% Adjust color limits for optimal performance 
for index=1:length(ah)
  axes(ah(index))
  shading faceted  
  caxis([-EMAX(index) +EMAX(index)])
end
kelicol
moveh(lx(2),-.015)
movev(lx,-.03)

% Create eigenvalue legends
movs=[repmat(-1.3,numf,1)...
      repmat(-1.3,numf,1)...
      repmat(-1.3,numf,1)];
for index=1:numf*length(m)
  axes(ah(index))
  lxb(index)=axes('position',lox(index,:));
  axis off
  ttl(index)=title(sprintf('%s = %9.6f','\lambda',...
			   Valo(index)));
  movev(ttl(index),movs(index))
end

fig2print(gcf,'portrait')
figdisp([],[],'-painters')

