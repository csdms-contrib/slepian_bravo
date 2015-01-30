function varargout=plm2slep(lmcosi,TH,L,phi,theta,omega,nosort,J)
% [falpha,V,N,MTAP,C]=PLM2SLEP(lmcosi,TH,L,phi,theta,omega,nosort,J)
%
% Finds the spherical harmonic expansion coefficients into a SINGLE-CAP,
% potentially rotated, Slepian basis of a function whose real spherical
% harmonic expansion coefficients are given.
%
% INPUT:
%
% lmcosi     Standard-type real spherical harmonic expansion coefficients
% TH         Radius of the concentration region (degrees) OR
%            'england', 'eurasia',  'namerica', 'australia', 'greenland'
%            'africa', 'samerica', 'amazon', 'orinoco', in which case 
%            you must have phi,theta,omega all equal to zero
%            OR: {'region' buf} where buf is the distance in degrees that
%            the region outline will be enlarged by BUFFERM
%            OR: [lon lat] an ordered list defining a closed curve [degrees]
% L          Bandwidth of the window [default: bandwidth of the data]
% phi        Longitude of the center of the tapers (degrees)
% theta      Colatitude of the center of the tapers (degrees)
% omega      Anticlockwise azimuthal rotation of the tapers (degrees)
% nosort     0 Will sort the output according to the global eigenvalue [default]
%            1 Will not sort thus the "block" sorting is preserved
% J          Number of largest eigenfunctions in which to expand [default: N]
%
% OUTPUT:
%
% falpha     The expansion coefficients of the function into the Slepian basis
% V          The eigenvalues of the Slepian functions in question
% N          The Shannon number
% MTAP       The orders of the Slepian functions in question, if preserved
% G          The matrix with the spherical harmonic expansion
%            coefficients of the Slepian functions used
%
% EXAMPLE:
%
% plm2slep('demo1') % Check that single-order functions correctly transform back
% plm2slep('demo2') % Slight variation on the above with multiple same-order ones
% plm2slep('demo3') % Moon without the South Pole Aitken Basin
%
% See also: PTOSLEP, GLMALPHA, GLMALPHAPTO, SLEP2PLM
%
% Last modified by fjsimons-at-alum.mit.edu, 06/26/2012
% With thanks to input from charig-at-princeton.edu, 06/26/2012

if ~isstr(lmcosi)
  % Supply defaults
  defval('TH',30)
  defval('L',18)
  defval('phi',0)
  defval('theta',0)
  defval('omega',0)
  defval('nosort',0)
  
  if lmcosi(1)~=0 || lmcosi(2)~=1
    error('Spherical harmonics must start from degree zero')
  end

  % If it is the standard North-Polar cap or a geographic region, it's easy
  if phi==0 && theta==0 && omega==0
    % Get the Slepian basis; definitely not block-sorted as for the rotated
    % versions this will make no sense at all anymore
    defval('J',round((L+1)^2*spharea(TH)))
    [G,V,EL,EM,N,GM2AL,MTAP,IMTAP]=glmalpha(TH,L,[],0,[],[],J);
  else
    % Need to get a complete GLMALPHA but for the rotated basis
    % Definitely, "single-order" has lost its meaning here, but the MTAP
    % will still identify what the order of the unrotated original was
    [G,V,EL,EM,N,GM2AL,MTAP,IMTAP]=glmalphapto(TH,L,phi,theta,omega);
  end

  if ~nosort
    % Sort by decreasing eigenvalue
    [V,vi]=sort(V,'descend');
    G=G(:,vi); if ~isnan(MTAP); MTAP=MTAP(vi); end
    % If you don't do this, the eigenfunctions are ordered in the way
    % that they correspond to single-orders back when, unrotated, they
    % belonged to a polar cap, and the eigenvalues are sorted within
    % these blocks. This is useful for, e.g. SPIE2009_1 a la SDSNEEUW.
  end

  % Get the mapping from LMCOSI into not-block-sorted GLMALPHA
  [~,~,~,~,~,~,~,~,~,ronm]=addmon(L);
  
  % Make sure that the requested L acts as truncation on lmcosi
  lmcosi=lmcosi(1:addmup(L),:);

  % Perform the expansion of the signal into the Slepian basis
  falpha=G'*lmcosi(2*size(lmcosi,1)+ronm(1:(L+1)^2));

  % Collect output
  varns={falpha,V,N,MTAP,G};
  varargout=varns(1:nargout);  
elseif strcmp(lmcosi,'demo1')
  % Change these parameters, keep alps below L-|m|+1
  TH=30; L=18; 
  % Pick whatever order comes to mind
  m=9; alps=7;
  m=(-1)^round(rand*2)*round(rand*L);
  alps=min(ceil(rand*L),L-abs(m)+1);
  
  warning('If the eigenvalues are tiny, sorting will be off and this misfires')
  
  disp(sprintf('\nChecking the %ith best function of order %i\n',alps,m))
  % Supply whatever rotation - but remember, GLMALPHAPTO is rerun for each
  phi=round(rand*360); theta=round(rand*180); omega=round(rand*360);
  phi=78; theta=78; omega=10;
  % Construct the alps first single-order rotated Slepian basis functions
  [lm,cosi]=ptoslep(phi,theta,omega,TH,L,m,alps);
  % Check the alps'th first
  [falpha,V,N,MTAP]=plm2slep([lm cosi(:,:,alps)],TH,L,phi,theta,omega);
  % Check that only one expansion coefficient is hit
  difer(falpha(indeks(find(MTAP==m),alps))-1)
elseif strcmp(lmcosi,'demo2')
  % Change these parameters, keep alps below L-|m|+1
  TH=30; L=18; 
  % Pick whatever order comes to mind
  m=9; alps=7;
  m=(-1)^round(rand*2)*round(rand*L);
  alps=min(round(rand*L),L-abs(m)+1);
  disp(sprintf('\nChecking the %ith best functions of order %i\n',alps,m))
  % Supply whatever rotation - but remember, GLMALPHAPTO is rerun for each
  phi=round(rand*360); theta=round(rand*180); omega=round(rand*360);
  phi=78; theta=78; omega=10;
  % Construct the alps first single-order rotated Slepian basis functions
  [lm,cosi]=ptoslep(phi,theta,omega,TH,L,m,alps);
  % Check all the alps first functins together
  [falpha,V,N,MTAP]=plm2slep([lm sum(cosi,3)],TH,L,phi,theta,omega);
  % Check that only the right expansion coefficients are hit
  difer(falpha(find(MTAP==m)-1))
elseif strcmp(lmcosi,'demo3')
  strin=lmcosi;
  % Do something for the Moon, quickly, from my own files
  lmcosi=fralmanac('GLTM2B','SHM'); 
  L=lmcosi(end,1); L=21
  % Figure windows
  clf
  [ah,ha]=krijetem(subnum(2,2));
  % Yup, it looks good
  axes(ah(1))
  [topo,b,c]=plotplm(lmcosi,[],[],4,1); delete([b c])
  % No let's take out the SPA basin, quickly spoken, that is
  TH=20;
  % I keep this small so I don't run into the pole - bear with me
  [lon2,lat2]=caploc([192 -61],TH,[]);
  hold on; d=plot(lon2,lat2,'k'); hold off

  % Color scale
  col1=getpos(ah(1));
  col1=[col1(1)-col1(3)/10 col1(2) col1(3)/15 col1(4)];
  [cb(1),xcb(1)]=addcb(col1,round(minmax(topo)),round(minmax(topo)),'jet');
  shrink(cb(1),1,1.55)
  
  % Let's make a big matrix with the Slepians for this closed coordinate
  % set; note that we are NOT using the fact that it is a circle so we
  % remain completely general. Though I fake the Shannon number ahead
  N1=round(spharea(TH,1)*(L+1)^2);
  % But this would have been about the same
  N2=round(spharea([lon2(:) lat2(:)])*(L+1)^2);
  % Compute the first so many Slepian functions of the anti-region
  anti=1;
  % How many will be needed? We take the Shannon number
  J=anti*(L+1)^2+(-1)^[anti]*N1;
  % This is an experimental parameter - take more and you delete less of
  % the structure but you'll get less mixture of the degrees hence the
  % other parts of the map remain more faithful
  
  % Make the filename to save them
  h=hash([lon2 lat2],'sha1');
  fname=fullfile(getenv('IFILES'),'GLMALPHA',...
		 sprintf('glmalpha-%s-%i-%i.mat',h,L,J));
  if anti==1
    % Update the file name to reflect the complimentarity of the region
    fname=sprintf('%s-1.mat',pref(fname)); 
  end
  
  % Get the Slepian eigenfunctions
  [G,V,EL,EM,N]=glmalpha([lon2 lat2],L,[],[],[],[],J,anti);
  
  % Now transform the regular spherical harmonic expansion to the Slepian ones
  % and then back to the spherical harmonics - but only for the
  % well-concentrated Slepian functions
  % Get the mapping from LMCOSI into not-block-sorted GLMALPHA
  [~,~,~,~,~,~,~,~,~,ronm]=addmon(L);
  stix=2*size(lmcosi,1)+ronm(1:(L+1)^2);
  
  % Perform the expansion of the signal into the Slepian basis
  lmcosiJ=lmcosi; lmcosiJ(:,3:4)=0;
  lmcosiJ(stix)=G*G'*lmcosi(stix);

  % Note that it is instructive to look at G*G' to see how the degrees
  % are being mixed to null the SPA structure

  % Yup, it looks good - you now have the spherical harmonic
  % representation of the lunar topography without the SPA basin, and you
  % can proceed with the regular analysis
  axes(ah(3))
  [topoJ,b,c]=plotplm(lmcosiJ,[],[],4,1); delete([b c])
  hold on; d=plot(lon2,lat2,'k'); hold off 

  % Color scale
  col2=getpos(ha(2));
  col2=[col2(1)-col2(3)/10 col2(2) col2(3)/15 col2(4)];
  [cb(2),xcb(2)]=addcb(col2,round(minmax(topoJ)),round(minmax(topoJ)),'jet');
  shrink(cb(2),1,1.55)

  % Take a look at the power spectrum
  axes(ah(2))
  [s,l]=plm2spec(lmcosi);
  p(1)=plot(l,decibel(s),'b-o'); axis tight
  axes(ah(4))
  [sJ,lJ]=plm2spec(lmcosiJ);
  p(2)=plot(lJ,decibel(sJ),'r-o'); axis tight
  
  % Cosmetics
  fig2print(gcf,'portrait')
  tits={sprintf('lunar topography to degree %i',L),...
	'power spectrum',...
	sprintf('lunar topography without SPA to degree %i',L)...
	,'power spectrum'};
  xes={'longitude','spherical harmonic degree',...
      'longitude','spherical harmonic degree'};
  yes={'latitude','spectral density (dB)',...
      'latitude','spectral density (dB)'};
  for in=1:length(ah)
    axes(ah(in))
    title(tits{in}); xlabel(xes{in}); ylabel(yes{in})
  end
  longticks(ah)
  set(ha(1:2),'yaxisl','r')
  set(ha(3:4),'yaxisl','r','xgrid','on','ygrid','on',...
	      'ylim',[-45 5],'xlim',[-5 L])
  shrink(ha(3:4),1,3/2)
  movev([ah(1:2) cb(1)],-.1)
  moveh([ha(1:2)],0.015)
  set(p(1),'MarkerS',2,'MarkerF','b')
  set(p(2),'MarkerS',2,'MarkerF','r')
  set(xcb,'string','')
  set(get(ha(1),'ylabel'),'rotation',-90)
  set(get(ha(2),'ylabel'),'rotation',-90)

  movev([ah cb],.15)
  
  hid=id;
  movev(hid,.15)
  
  figna=figdisp([],strin,[],1);
  system(sprintf('degs %s.eps',figna));
  system(sprintf('epstopdf %s.eps',figna));

  % After this you would simply be doing the PLM2ROT analysis of the
  % first few coefficients

  %   The coordinates for the SPA deselection we chose to be a center of
  %   -55 deg. lat and 10 deg. lon for a map centered on the farside with a
  %   radius of 1000 km. There is flexibility though in this choice. The
  %   SPA deselection should happen on 'global' fits as well as farside
  %   hemisphere-only fits. The third fit type is nearside hemisphere
  %   only. 
end
