function varargout=...
    xyz2slep(fthph,theta,phi,TH,L,phi0,theta0,omega,J,sord,Glma,V,N,EL,EM)
% [falpha,N,V,Glma,EL,EM,lmcosi]=...
%    XYZ2SLEP(fthph,theta,phi,TH,L,phi0,theta0,omega,J,sord,Glma,V,N,EL,EM)
%
% Slepian expansion of irregularly sampled but closely collocated points
% with respect to a Slepian basis that is bandlimited and spatially
% concentrated to a (rotated) single cap of some description.
%
% INPUT:
%
% fthph       function values defined on the array theta,phi
% theta       colatitude vector (0 <= theta <= pi) [radians]
% phi         longitude vector (0 <= theta <= 2*pi) [radians]
% TH          the radius of the single cap Slepian function
% L           the bandlimit of the Slepian function, or the passband
% phi0        the longitude of the center of the cap [degrees] 
% theta0      the latitude of the center of the cap [degrees]
% omega       the rotation around the cap center [degrees]
% J           the truncation number in the expansion [default: N]
% sord        1 Single cap of diameter 2TH [default]
%             3 Equatorial belt of width 2TH
% Glma        The spectral eigenfunctions in case you already have them
% V           The spectral eigenvalues in case you already have them
% N           The Shannon number in case you already have it
% EL          The degrees in question if you already have them
% EM          The orders in question if you already have them
%
% OUTPUT:
%
% falpha      the J expansion coefficients in the new basis
% N           the Shannon number
% V           the eigenvalue sequence
% Glma        the spectral eigenfunctions, for use in, e.g. PLOTSLEP
% V           the spectral eigenvalues
% N           the Shannon number
% EL          the degrees in question
% EM          the orders in question
% lmcosi      the spherical harmonic expansion from XYZ2PLM
%
% EXAMPLE:
%
% xyz2slep('demo1',N) with optional N number of samples
% xyz2slep('demo2',N) with optional N number of samples
% xyz2slep('demo3',N) with optional N number of samples
% 
% SEE ALSO:
%
% PLM2SLEP, XYZ2SPL, XYZ2PLM
%
% Last modified by fjsimons-at-alum.mit.edu, 01/17/2010

% Later: modify to do double cap nonrotated, and compliment

% Supply some default values for the Slepian basis
if ~isstr(fthph)
  defval('TH',15);
  defval('L',36);
  defval('phi0',15);
  defval('theta0',70);
  defval('omega',0);
  defval('sord',1);
  
  % Supply some default values for the function values
  defval('Nd',300);
  [lon,lat]=randpatch(Nd,TH,theta0,phi0);
  defval('theta',pi/2-lat*pi/180)
  defval('phi',lon*pi/180);
  
  % Construct the spatial eigenfunctions evaluated on this irregular patch
  if ~all(size(theta(:))==size(phi(:)))
    error('Input arrays must have the same dimensions for irregular grids')
  end
  % Compute the Shannon number for the relevant situation
  N=round((L+1)^2*spharea(TH,sord));
  % And equate the truncation level to the Shannon number unless
  % otherwise indicated  
  defval('J',round(N))

  if sord==1
      if exist('Glma')~=1 ||  exist('V')~=1 ||  exist('N')~=1 ...
	|| exist('EL')~=1 ||  exist('EM')~=1
	if phi0~=0 || theta0~=0 || omega~=0
	  % (Rotated) single cap
	  [Gar,V,N,J,phi0,theta0,omega,theta,phi,TH,L,Glma,EL,EM]=...
	      galphapto(TH,L,phi0,theta0,omega,theta,phi,J,1);
	else
	  % (Non-Rotated) single cap
	  [Gar,V,EM,GK,VK,NA,N,theta,phi,Glma,EL]=...
	      galpha(TH,L,1,theta,phi,'global',0,0,0,J,1);
	end
      else
	if phi0~=0 || theta0~=0 || omega~=0
	  Gar=galphapto(TH,L,phi0,theta0,omega,theta,phi,J,1,Glma,V,N,EL,EM);
	else
	  Gar=galpha(TH,L,1,theta,phi,'global',0,0,0,J,1,Glma,V,N,EL,EM);
	end
	end
  elseif sord==3
    meths='new';
    switch meths
     case 'old'
      % Polar double cap of 90-TH 
      [Gar,V]=galpha(90-TH,L,2,theta,phi,'global',0,0,[],[],1);
      % Belt of width 2TH
      [a,b]=sort(V);
      V=1-sort(V);
      % Don't just flipud as the eigenvalues are very close
      Gar=Gar(b,:);
      Gar=Gar(1:J,:);
      V=V(1:J);
     otherwise
       % WHICH SHOULD BE EQUIVALENT TO THE NEW
       [Gar2,V2]=galpha(TH,L,3,theta,phi,'global',0,0,[],J,1);
    end
  end
  Nd=length(theta(:));
  
  if J<Nd
    % Overdetermined problem, construct the least-squares solution using pinv
    falpha=pinv(Gar')*fthph;
    % Should be well-conditioned, as the eigenvalues are already, i.e. it
    % should be identical to inv(Gar*Gar')*Gar, which it is
    % truncated. Quality should be roughly according to how close G'G
    % really is to V.
  else
    error('Think about this')
  end

  if nargout>6
    % Test this... for fun
    disp('I am doing it also in spherical harmonics!')
    lmcosi=xyz2plm(fthph,L,[],90-theta*180/pi,phi*180/pi);
  end
  
  % Initialize variables that could have ended empty
  defval('Glma',NaN);
  defval('EL',NaN);
  defval('EM',NaN);
  defval('lmcosi',NaN);
  
  % Prepare output
  vars={falpha,N,V,Glma,EL,EM,lmcosi};
  varargout=vars(1:nargout);
  
elseif strcmp(fthph,'demo1')
  defval('theta',[])
  defval('Nd',theta);
  defval('Nd',300)
  theta=[];
  % What if the input is a single sampled Slepian function?
  TH=15;
  phi0=15;
  theta0=70;
  L=36;
  [lon,lat]=randpatch(Nd,TH,phi0,theta0);  
  defval('theta',pi/2-lat*pi/180)
  defval('phi',lon*pi/180);
  % Generate the data at the sample points
  [Gar,V,N,J]=galphapto(TH,L,phi0,theta0,omega,theta,phi,[],1);
  alfa=ceil(rand*J);
  % Pick the alfa'th basis function as a first example
  fthph=Gar(alfa,:)';
  
  % Add some noise?
  nadd=rand(size(fthph))*std(fthph)/20;
  SN=rms(fthph)/rms(nadd);
  fthph=fthph+nadd;

  % Now attempt to invert and find the single coefficient back
  falpha=xyz2slep(fthph,theta,phi,TH,L,phi0,theta0,[],J);
  
  disp(sprintf('error in alpha %8.3f',max(falpha)-1))
  disp(sprintf('error in not alpha %8.3f',mean(abs(skip(falpha,alfa)))))
  
  % Report on the performance
  testreport(falpha,Gar,fthph)
  
  % Plot the "original" data from which sampled on an appropriate grid
  [Gar,V,N,J,phi0,theta0,omega,thetag,phig]=galphapto(TH,L,phi0,theta0);

  % Define the plot boundaries
  c11cmn=[phi0 90-theta0 phi0 90-theta0]+2*[-TH TH TH -TH];

  clf
  [ah,ha]=krijetem(subnum(2,2));
  axes(ah(1))
  data=reshape(Gar(alfa,:)',length(thetag),length(phig));
  colaxis=max(abs(data(:)))*[-1 1];
  data1=data;
  
  % Cut out the smallest values for ease of visualization
  imagefnan(c11cmn(1:2),c11cmn(3:4),setnans(data1),[],colaxis)
  ploco(c11cmn,theta0,phi0,TH,'original',data)

  axes(ah(2))
  e=plot(lon-360*[lon>180],lat,'kv','MarkerF','w');
  ploco(c11cmn,theta0,phi0,TH)
  title(sprintf('N = %i, N/S = %3.1f%s',length(fthph),100/SN,'%s'))

  axes(ah(3))
  data2=reshape(Gar'*falpha,length(thetag),length(phig));
  data3=data2;
  
  % Cut out the smallest values for ease of visualization
  imagefnan(c11cmn(1:2),c11cmn(3:4),setnans(data3),[],colaxis)
  ploco(c11cmn,theta0,phi0,TH,'reconstruction',data2,data)

  axes(ah(4))
  data4=data-data2;
  imagefnan(c11cmn(1:2),c11cmn(3:4),data4,[],colaxis,[],colaxis)
  ploco(c11cmn,theta0,phi0,TH,'error',data4,data)
  
  [cb,xl]=addcb('vert',round(colaxis),round(colaxis),'kelicol',1,1);
  set(xl,'string','reconstruction error')
  set(cb,'YAxisL','r')

  fig2print(gcf,'portrait')
elseif strcmp(fthph,'demo2')
  defval('theta',[])
  defval('Nd',theta);
  defval('Nd',300)
  theta=[];
  % What if the input is a sampled set of Slepian functions?
  TH=15;
  phi0=15;
  theta0=70;
  L=36;
  [lon,lat]=randpatch(Nd,TH,phi0,theta0);  
  defval('theta',pi/2-lat*pi/180)
  defval('phi',lon*pi/180);
  % Generate the data
  [Gar,V,N,J]=galphapto(TH,L,phi0,theta0,[],theta,phi,[],1);

  % Pick a random combination of basis functions 
  salpha=randn(J,1);
  fthph=Gar'*salpha;
  
  % Add some noise?
  nadd=rand(size(fthph))*std(fthph)/20;
  SN=rms(fthph)/rms(nadd);
  fthph=fthph+nadd;

  % Now attempt to invert and find the coefficient back
  [falpha,N]=xyz2slep(fthph,theta,phi,TH,L,phi0,theta0,[],J);

  disp(sprintf('rmse in alpha %3.1f%s',...
	       rms(falpha-salpha)/rms(salpha)*100,'%'))
  
  % Report on the performance
  testreport(falpha,Gar,fthph)

  % Plot the "original" data from which sampled on an appropriate grid
  [Gar,V,N,J,phi0,theta0,omega,thetag,phig]=galphapto(TH,L,phi0,theta0);

  % Define the plot boundaries
  c11cmn=[phi0 90-theta0 phi0 90-theta0]+2*[-TH TH TH -TH];
  
  clf
  [ah,ha]=krijetem(subnum(2,2));
  axes(ah(1))
  data=reshape(Gar'*salpha,length(thetag),length(phig));
  colaxis=max(abs(data(:)))*[-1 1];
  data1=data;
  
  % Cut out the smallest values for ease of visualization
  imagefnan(c11cmn(1:2),c11cmn(3:4),setnans(data1),[],colaxis)
  ploco(c11cmn,theta0,phi0,TH,'original',data)

  axes(ah(2))
  e=plot(lon-360*[lon>180],lat,'kv','MarkerF','w');
  ploco(c11cmn,theta0,phi0,TH)
  title(sprintf('N = %i, N/S = %3.1f%s',length(fthph),100/SN,'%s'))
  axes(ah(3))
  data2=reshape(Gar'*falpha,length(thetag),length(phig));
  data3=data2;
  
  % Cut out the smallest values for ease of visualization
  imagefnan(c11cmn(1:2),c11cmn(3:4),setnans(data3),[],colaxis)
  ploco(c11cmn,theta0,phi0,TH,'reconstruction',data2,data)

  axes(ah(4))
  data4=data-data2;
  imagefnan(c11cmn(1:2),c11cmn(3:4),data4,[],colaxis)
  ploco(c11cmn,theta0,phi0,TH,'error',data4,data)

  [cb,xl]=addcb('vert',round(colaxis),round(colaxis),'kelicol',10,1);
  set(xl,'string','reconstruction error')
  set(cb,'YAxisL','r')

  fig2print(gcf,'portrait')
elseif strcmp(fthph,'demo3')
  defval('theta',[])
  defval('Nd',theta);
  defval('Nd',300)
  theta=[];

  TH=15;
  L=36;
  phi0=15;
  theta0=70;
  [lon,lat]=randpatch(Nd,TH,phi0,theta0);  
  defval('theta',pi/2-lat*pi/180)
  defval('phi',lon*pi/180);

  % Generate the data
  lmcosi=load(fullfile(getenv('IFILES'),...
		       'EARTHMODELS','POMME-4','pomme-4.2s-nosecular.cof'));
  % Restrict to degree L
  lmcosi=lmcosi(1:addmup(L)-addmup(lmcosi(1)-1),:);

  % Convert to radial-component magnetic field on the reference surface
  lmcosip=plm2mag(lmcosi);
  
  % Expand at the random set of points  
  fthph=plm2xyz(lmcosip,lat,lon);
  
  % Add some noise?
  nadd=rand(size(fthph))*std(fthph)/2;
  SN=rms(fthph)/rms(nadd);
  fthph=fthph+nadd;

  % Find the appropriate Shannon number for this patch
  J=round((L+1)^2*spharea(TH,1));
  
  % Now attempt to make the Slepian expansion of the input
  [salpha,V1]=plm2slep([0 0 0 0 ; lmcosip],TH,L,phi0,theta0);

  % Now attempt to invert and find the coefficient back
  [falpha,N,V2]=xyz2slep(fthph,theta,phi,TH,L,phi0,theta0,[],J);
  
  difer(V1(1:J)-V2)
  difer(J-round(N))
  
  rms(salpha(1:J)-falpha)/rms(salpha(1:J))*100
  
  % Find the Slepian functions at the sampled points
  Gar=galphapto(TH,L,phi0,theta0,[],theta,phi,[],1);

  % Report on the performance
  testreport(falpha,Gar,fthph)

  % Define the plot boundaries
  c11cmn=[phi0 90-theta0 phi0 90-theta0]+2*[-TH TH TH -TH];

  clf
  [ah,ha]=krijetem(subnum(2,2));

  % Plot the original data on a complete grid
  axes(ah(1))
  data=plm2xyz(lmcosip,1,c11cmn);
  colaxis=max(abs(data(:)))*[-1 1];
  data1=data;
  
  % Cut out the smallest values for ease of visualization
  imagefnan(c11cmn(1:2),c11cmn(3:4),setnans(data1),[],colaxis)
  ploco(c11cmn,theta0,phi0,TH,'original',data)

  axes(ah(2))
  e=plot(lon-360*[lon>180],lat,'kv','MarkerF','w');
  ploco(c11cmn,theta0,phi0,TH)
  title(sprintf('N = %i, N/S = %3.1f%s',length(fthph),100/SN,'%s'))

  % Plot the inversion results on a complete grid
  axes(ah(3))
  [Gar,V,N,J,phi0,theta0,omega,thetag,phig]=galphapto(TH,L,phi0,theta0);
  data2=reshape(Gar'*falpha,length(thetag),length(phig));
  data3=data2;
  
  % Cut out the smallest values for ease of visualization
  imagefnan(c11cmn(1:2),c11cmn(3:4),setnans(data3),[],colaxis)
  ploco(c11cmn,theta0,phi0,TH,'reconstruction',data2,data)

  axes(ah(4))
  data4=data-data2;
  imagefnan(c11cmn(1:2),c11cmn(3:4),data4,[],colaxis)
  ploco(c11cmn,theta0,phi0,TH,'error',data4,data)

  [cb,xl]=addcb('vert',round(colaxis),round(colaxis),'kelicol',20000,1);
  set(xl,'string','reconstruction error')
  set(cb,'YAxisL','r')

  fig2print(gcf,'portrait')
end

% Some auxiliary functions useful for the demos above
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ploco(c11cmn,theta0,phi0,TH,strng,data,dataref)
% Plot the continents - note the Greenwich trick
yesorno=360*[c11cmn(1)<0 0];
[ax,f,lola1]=plotcont(c11cmn(1:2)+yesorno,...
		      c11cmn(3:4)+yesorno,[],-360);
[ax,f,lola2]=plotcont([0 c11cmn(2)],c11cmn(3:4));
axis(c11cmn([1 3 4 2]))
% Plot the circle of concentration
hold on
[lon2,lat2]=caploc([phi0 90-theta0],TH);
f=plot(lon2-360*[lon2>180],lat2,'k-');
% Plot the center
d=plot(phi0,90-theta0,'o','MarkerF','w','MarkerE','k');
if nargin>=6
  % Figure out the data INSIDE versus the data OUTSIDE
  latd=90-[theta0-2*TH:1:theta0+2*TH];
  lond=[phi0-2*TH:1:phi0+2*TH];
  if ~all([length(latd)==size(data,1) length(lond)==size(data,2)])
    error('Wrong assumption about the data size')
  end
  % Distance from center
  [LOND,LATD]=meshgrid(lond,latd);
  [gcdkm,delta]=grcdist([phi0 90-theta0],[LOND(:) LATD(:)]);
  rrms=sqrt(mean(data(delta<=TH).^2));
  if nargin==7
    refrms=sqrt(mean(dataref(delta<=TH).^2));
    title(sprintf('%s, region-rms = %8.3f or %3.1f%s',...
		  strng,rrms,rrms/refrms*100,'%'))
  else
    title(sprintf('%s, region-rms = %8.3f',...
		  strng,rrms))
  end
end
set(gca,'xtick',phi0+[-2*TH -TH 0 TH 2*TH],...
	'ytick',90-theta0+[-2*TH -TH 0 TH 2*TH])
grid on
hold off
drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data=setnans(data)
data(abs(data)<max(abs(data(:)))/1000)=NaN;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function testreport(falpha,Gar,fthph)
% How are the data being reconstructed at the sampled points?
fthphest=Gar'*falpha;

% Comment on the rms of the misfit
rmsf=sqrt(mean((fthph).^2));
rmse=sqrt(mean((fthph-fthphest).^2));
disp(sprintf('\nrms  of the sampling points = %8.3f',...
	     rmsf))
disp(sprintf('rmse at the sampling points = %8.3f or %3.1f%s\n',...
	     rmse,rmse/rmsf*100,'%'))
  
