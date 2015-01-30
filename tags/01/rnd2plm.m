function varargout=rnd2plm(lmcosi,meth)
% lmcosip=RND2PLM(lmcosi,meth)
%
% Randomizes a spherical-harmonic coefficient matrix
%
% INPUT:
%
% lmcosi    Input matrix listing l, m, Ccos and Csin
% meth      'kevin' a la Kevin (faster for low degrees)
%           'frederik' a la Frederik (faster for high degrees)
%
% OUTPUT:
%
% lmcosip   Output matrix listing l, m, Ccos and Csin
%
% SEE ALSO: 
%
% PLM2RND, COV2PLM
% 
% REQUIRES:
%
% RandOrthMat, from
% http://www.mathworks.com/matlabcentral/fileexchange/11783-randorthmat
%
% EXAMPLE: 
%
% RND2PLM('demo1') produces realizations of fields with the same
% statistics as lunar topography
%
% Written by Kevin Lewis, 02/21/2010
% Last modified by fjsimons-at-alum.mit.edu, 07/06/2012

if ~isstr(lmcosi)
  % Supply the default, i.e. a low-resolution terrestrial gravity field
  defval('lmcosi',...
	 kindeks(rindeks(fralmanac('EGM96','SHM'),1:addmup(255)-3),1:4))
  defval('meth','kevin')
  
  switch meth 
   case 'kevin'
    tic
    l=lmcosi(:,1);
    m=lmcosi(:,2);
    co=lmcosi(:,3);
    si=lmcosi(:,4);
    for i=min(l):max(l)
      ind=find(l==i);
      lcoeffs=RandOrthMat(2*i+1)*[flipud(si(ind(2:end))) ; co(ind)];
      si2(ind(2:end))=flipud(lcoeffs(1:i));
      co2(ind)=lcoeffs(i+1:end);
    end
    lmcosip=[l m co2(:) si2(:)];
    toc
   case 'frederik'
    tic
    lmcosip=[lmcosi(:,1:2) zeros(length(lmcosi),2)];
    % Loop over the degrees
    for l=lmcosi(1,1):lmcosi(end,1)
      % Extract the COSINE coefficients and their indices at increasing degree
      [Ccos,b,e]=shcos(lmcosi,l);
      % Extract the SINE coefficients and return them at decreasing degree
      Csin=lmcosi(e:-1:b+1,4);
      % Randomize by orthogonal multiplication
      lcoeffs=RandOrthMat(2*l+1)*[Csin ; Ccos];
      % Assign the SINE coefficients (m>0) to the right spots in the larger matrix
      % Note the trick with the sum to convert the empty to a zero if needed
      lmcosip(e:-1:b+1,4)=sum(lcoeffs(1:l),2);
      % Assign the COSINE (m<=0) coefficients to the right spots
      lmcosip(b:e,3)=lcoeffs(l+1:end);
    end
    toc
  end

  % Prepare output
  varns={lmcosip};
  varargout=varns(1:nargout);
  
elseif strcmp(lmcosi,'demo1')
  clf
  % Capture the demo string
  strin=lmcosi;
  % Load lunar topography
  lmcosi=fralmanac('GLTM2B','SHM'); 
  L=lmcosi(end,1); 
  [ah,ha]=krijetem(subnum(3,2));
  
  axes(ah(1))
  [topo,b,c]=plotplm(lmcosi,[],[],4,1); delete([b c])
  % Color scale
  col1=getpos(ah(1));
  col1=[col1(1)-col1(3)/10 col1(2) col1(3)/15 col1(4)];
  [cb(1),xcb(1)]=addcb(col1,round(minmax(topo)),round(minmax(topo)),'jet');
  shrink(cb(1),1,1.55)
  
  % Calculate the global power spectral density
  [SWSl,l]=plm2spec(lmcosi);

  axes(ah(2))
  p(1)=plot(l,decibel(SWSl),'b-o'); axis tight
  
  % Prepare for the localization
  TH=20;
  Lw=10;
  phi0=200;
  th0=90;
  
  axes(ha(3))
  [G,V,EL,EM,N]=glmalphapto(TH,Lw,phi0,th0,0);
  [V,i]=sort(V,'descend'); G=G(:,i); 
  wot=1;
  % Check the normalization of GLM2LMCOSI
  lmcosit=glm2lmcosi(G,wot);
  [taper,b,c]=plotplm(lmcosit,[],[],4,1); delete([b c])
  
  [st,lt]=plm2spec(lmcosit);
  axes(ah(6))
  p(2)=plot(lt,decibel(st),'b-o'); axis tight

  % No point in generating more global examples as the global power
  % spectrum is always identical 
  for inde=1:25
    % Generate random topographies just like the moon
    lmcosip=rnd2plm(lmcosi);
    topop=plm2xyz(lmcosip,1);
  
    % Plot examples in the space domain
    axes(ah(3))
    [~,b,c]=plotplm(topop,[],[],4,1); delete([b c])
    
    [lon2,lat2]=caploc([phi0 90-th0],TH,[]);
    hold on; d=plot(lon2,lat2,'k'); hold off
    
    col2=getpos(ha(2));
    col2=[col2(1)-col2(3)/10 col2(2) col2(3)/15 col2(4)];
    [cb(2),xcb(2)]=addcb(col2,round(minmax(topop)),round(minmax(topop)),'jet');
    shrink(cb(2),1,1.55)
    
    Sal=0;
    % Calculate the LOCAL power spectral density
    for onde=1:min(round(3*N),size(G,2))
      % To keep with Dahlen & Simons, need the factor
      % See RB VIII, p 126 
      taper=plm2xyz(glm2lmcosi(G,onde),1)*sqrt(4*pi);
      % Note that the XYZ2PLM uses 4pi normalized harmonics
      lmcosipt=xyz2plm(topop.*taper,L);
      lmcosipt(:,3:4)=lmcosipt(:,3:4)/sqrt(4*pi);
      Sal=Sal+V(onde)*plm2spec(lmcosipt);
    end
    SMTl=Sal/sum(V(1:onde));
    axes(ah(4))
    p(2+inde)=plot(0:L,SMTl,'b-o'); axis tight
    hold on
    
    drawnow
  end
  
  % So now the claim is that the AVERAGE of the SMTl is given by the
  % convolution of the multitaper coupling kernel with the true global
  % spectrum and that the VARIANCE of the SMTl is given by the multitaper
  % variance (which is in function of the truth).

  % Now calculate the expected value of the multitaper estimates in
  % function of the truth. Took 80 because 70 has a file problem
  Mllp=mcouplings(Lw,80);
  Mllp=Mllp(1:L+1,1:L+1);
  
  ESMTl=Mllp*SWSl;

  % Calculate the variance-use the global variance??
  varSl=mtvar(SWSl,l,Lw,TH,1);
  [pp,ex,ey]=errorxy(l,ESMTl,[],2*sqrt(varSl));
  hold on
  plot(l,ESMTl,'k-','LineW',2)
  
  keyboard
  
  % Figure cosmetics
  fig2print(gcf,'tall')
  tits={sprintf('lunar topography to degree %i',L),...
	'power spectrum',...
	sprintf('randomized lunar topography %i',L)...
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
  set(ha(1:3),'yaxisl','r')
  set(ha(4:6),'yaxisl','r','xgrid','on','ygrid','on',...
	      'ylim',[-45 5],'xlim',[-5 L])
  shrink(ha(3:4),1,3/2)
  movev([ah(1:2) cb(1)],-.1)
  moveh([ha(1:2)],0.015)
  set(p(1:2),'MarkerS',2,'MarkerF','b')
  set(p(3:end),'MarkerS',2,'MarkerF','r')
  set(xcb,'string','')
  set(get(ha(1),'ylabel'),'rotation',-90)
  set(get(ha(2),'ylabel'),'rotation',-90)
  movev([ah cb],.15)
  
  figna=figdisp([],sprintf('%s_%i',strin,inde),[],1);
  system(sprintf('degs %s.eps',figna));
  system(sprintf('epstopdf %s.eps',figna));
  
  keyboard
end
