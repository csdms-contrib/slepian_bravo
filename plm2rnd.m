function varargout=plm2rnd(L,bta,meth,norma)
% [lmcosi,bta,bto,sdl,el]=PLM2RND(L,BTA,meth,norma)
%
% Makes a random field up to degree L with a power spectral slope scaled
% to exactly l^bta. 
% 
% INPUT:
%
% L         Maximum degree of the expansion
% bta       Spectral slope (default: geopotential Earth: -4.0361)
% meth      1 Start from Gaussian random variates
%           2 Start from uniform random variates
% norma     The normalization used in PLM2SPEC (default: 2)
%           1 *(l+1) 
%           2 /(2l+1)
%           3 none
%
% OUTPUT:
%
% lmcosi    Matrix with [l m Ccos Csin] spherical harmonics
% bta       Input beta value, again
% bto       Actual spectral slope of this realization
% sdl       Actual spectrum of this realization
% el        Degrees of the spectrum
%
% SEE ALSO:
%
% RND2PLM, COV2PLM
%
% EXAMPLE:
%
% plm2rnd('demo1') % Verifies input and output
% plm2rnd('demo2') % Compares EGM96 to random structures
%
% Last modified by fjsimons-at-alum.mit.edu, 03/13/2013

defval('bta',-4.0361)
defval('L',100)
defval('meth',1)
defval('norma',2)

if ~isstr(L)
  [m,l,mzero]=addmon(L);
  switch meth
   case 1
    disp('Using Gaussian random')
    % One realization for every coefficient
    c=randn(addmup(L),2);
    % Make sure the sine coefficients are zero
    c(mzero,2)=0; 
    % Collect in a proper format
    lmcosi=[l m c];
   case 2
    disp('Using uniform random')
    c=2*rand(addmup(L),2)-1;
    c(mzero,2)=0;
    lmcosi=[l m c]; 
   otherwise
    error('Specify valid method')
  end
  
  % Compute the spectrum defined by plm2spec
  [sdl,el,bto]=plm2spec(lmcosi,norma);
  % Repeat observed spectrum for all m
  srep=addmin(sdl);
  % Then reassign it into the coefficient matrix after scaling
  lmcosi(:,3)=lmcosi(:,3)./sqrt(srep).*l.^(bta/2);
  lmcosi(:,4)=lmcosi(:,4)./sqrt(srep).*l.^(bta/2);
  % Fix the mean to be zero on average
  lmcosi(1,3)=2*round(rand)-1;
  lmcosi(1,4)=0;
  
  % And then see what you get
  if nargout>=3
    [sdl,el,bto]=plm2spec(lmcosi,norma);
    disp(sprintf('Best-fitting beta %5.2f',bto))
  end
  varns={lmcosi,bta,bto,sdl,el};
  varargout=varns(1:nargout);
else
  switch L
   case 'demo1'
    L=round(rand*180);
    [lmcosi,bta1]=plm2rnd(L,[]);
    [psd,l,bta2,lfit,logy,logpm]=plm2spec(lmcosi,2);
    disp(sprintf('L= %i ; slope in= %8.3f ; slope out= %8.3f',...
		 L,bta1,bta2))
   case 'demo2'    
    egm=fralmanac('EGM96','SHM');
    ah=krijetem(subnum(3,2));
    axes(ah(1))
    plotplm(egm(4:end,:),[],[],1)
    axes(ah(2))
    dbt=plotplm(egm,[],[],3);

    syn1=plm2rnd(max(egm(:,1)),bta);
    syn2=plm2rnd(max(egm(:,1)),bta);
    
    axes(ah(3))
    plotplm(syn1(7:end,:),[],[],1)
    
    axes(ah(4))
    dbt1=plotplm(syn1(4:end,:),[],[],3);

    axes(ah(5))
    plotplm(syn2(7:end,:),[],[],1)
    axes(ah(6))
    dbt2=plotplm(syn1(4:end,:),[],[],3);
    figdisp([],'demo2')
  end
end

% Some checks and balances
jj=1;
if jj==0
  % This would be like xyz2rnd or rnd2xyz
  fax=3;
  sdlav=0;
  N=10;
  % All of this should be better on Fibonacci grids
  for in=1:N
    r=randn(fax*38,fax*74);
    r(1,:)=r(1,randi(size(r,2)));
    r(end,:)=r(end,randi(size(r,2)));
    r(:,end)=r(:,1);
    % Average of the signal squared over the area in the xyz2plm
    % normalization is easily obtained by (not varying much with j)

    % Attempt to relate (1/4/pi)*int(signal^2)dOmega
    % Integral version of the power (energy per area)
    % Only want the first one since I'm not looking for higher moments
    j=20; lmcosi=xyz2plm(r.^2,j); a=lmcosi(1,3);
    % Discrete version of the power (energy per area)
    b=var(r(:)); 
    % Total power via Parseval - try high-end degrees, use integration
    j=200;
    warning off
    lmcosi=xyz2plm(r,j,'gl');
    warning on
    c=sum(sum(lmcosi(:,3:4).^2)); 
    [sdl,l]=plm2spec(lmcosi);
    % Total power (energy per area)
    d=sum(sdl.*(2*l+1));
    [a b c d]
    
    % Keep track of the average
    sdlav=sdlav+sdl;
    hold on
    plot(l,sdl)
    drawnow
  end
  sdlav=sdlav/N;
  plot(l,sdlav,'y-','linew',2)
  xlabel('spherical harmonic degree')
  ylabel('power spectral density')
  % Should close this loop, what should this look like??
  plot([l(1) l(end)],repmat(1/(j+1)^2,1,2),'m-')
  % But the constant part is here
  j=80;
  plot([l(1) l(j+1)],repmat(mean(sdlav(1:j+1)),1,2),'r-','Linew',2)
  longticks(gca)
  hold off
end
