function varargout=xyz2spl(flatlon,lat,lon,latp,lonp,method,pars)
% flatlonp=xyz2spl(flatlon,lat,lon,latp,lonp,method,pars)
%
% Spline interpolation of irregularly sampled spherical functions
%
% INPUT:
%
% fthph         Function values at given lat and lon
% lat           Given latitudes, in degrees
% lon           Given longitudes, in degrees
% latp          Desired latitudes, in degrees
% lonp          Desired longitudes, in degrees
% method        'abelpoisson' [default]
% pars          The parameter, e.g. h for Abel-Poisson [default: 0.9]
%
% OUTPUT:
%
% flatlonp      Function values at desired lat and lon
%
% EXAMPLES:
%
% XYZ2SPL('demo1')
% XYZ2SPL('demo2')
%
% Should also do an example with the competing method
%
% Last modified by fjsimons-at-alum.mit.edu, 05/16/2008

defval('flatlon',0)

if ~isstr(flatlon)
  % Specify what kind of spline and the parameters involved
  defval('lat',0)
  defval('lon',0)
  defval('method','abelpoisson')
  defval('pars',0.90)
  if sum(flatlon.*lat.*lon)==0
    % Bogus example
    [lon,lat]=randsphere(100);
    flatlon=rand(100,1);
    lonp=lon; latp=lat;
  end
  
  % Make sure they are all column vectors
  flatlon=flatlon(:);
  
  % So, what we need is the Gram matrix with the relevant spline kernel
  % Convert to colatitude and longitude in radians
  [TH,THp]=meshgrid((90-lat(:))*pi/180);
  [PH,PHp]=meshgrid(lon(:)*pi/180);
  
  % Make the matrix with the epicentral distances of all the combinations
  cosDelta=cos(TH).*cos(THp)+sin(TH).*sin(THp).*cos(PH-PHp);

  % And of the desired locations with the old locations, so no longer square
  [TH,THPp]=meshgrid((90-latp(:))*pi/180,(90-lat(:))*pi/180);
  [PH,PHPp]=meshgrid(lonp(:)*pi/180,lon(:)*pi/180);
  cosDeltaP=cos(TH).*cos(THPp)+sin(TH).*sin(THPp).*cos(PH-PHPp);
  
  switch method
   case 'abelpoisson'
    h=pars;
    % Calculate the Abel-Poisson kernel for these pairwise points
    K=(1-h^2)/4/pi./(1+h^2-2*cosDelta*h).^(3/2);
    KP=(1-h^2)/4/pi./(1+h^2-2*cosDeltaP*h).^(3/2);
  end
   
  % Now find the expansion coefficients for the spline
  ak=pinv(K)*flatlon;
  
  % Return the function interpolated at the desired longitudes etc
  % Interpolation to the same points would be identical
  flatlonp=[ak(:)'*KP]';
  
  % Provide output
  vars={flatlonp};
  varargout=vars(1:nargout);
elseif strcmp(flatlon,'demo1')
  % Inspired from PLM2XYZ('demo2')
  % We make up a function, sample randomly from it, and see how it
  % approximates the underlying function while fitting the sampled points
  % exactly 
  lmax=10; L=10;
  [m,l,mzero]=addmon(lmax);
  c=randn(addmup(lmax),2).*([l l].^(-1));
  c(1)=3; c(mzero,2)=0; lmcosi=[l m c];
  [r,lon,lat]=plm2xyz(lmcosi,180/sqrt(L*(L+1)));
  tol=length(lon)*length(lat);
  defval('degres',0.4)
  fra=degres;
  unform=2;
  [LON,LAT]=meshgrid(lon,lat);
  % Really uniform on the sphere
  [lonr,latr]=randsphere(100);%ceil(fra*tol));
  indo=sub2ind(size(r),ceil(scale(latr,[1 length(lat)])),...
		 ceil(scale(lonr,[1 length(lon)])));
  latr=LAT(indo);
  
  % Decide on the width of an Abel-Poisson kernel by choosing h
  h=rand(1);

  % Now use the function here for r(indo) and reconstruct it on a whole
  % fine grid... from lon(INDO) back to lon
  flatlonp=xyz2spl(r(indo),LAT(indo),LON(indo),LAT(:),LON(:),...
		   [],h);
  rp=reshape(flatlonp,length(lat),length(lon));
  
  % Predecide on an y-axis range
  ylix=1.1*minmax([r(:) rp(:)]);
  
  % Now we verify that the ones that we provided are fit EXACTLY
  difer(flatlonp(indo)-r(indo))
  
  % And we judge visually how we've smoothed around the unknown locations
  clf
  [ah,ha]=krijetem(subnum(2,2));

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  axes(ah(1))
  % Plot the original function in Cartesian projection
  plotplm(r,lon*pi/180,lat*pi/180,4)
  % Plot the sample locations
  hold on
  % [x,y]=mollweide(LON(indo)*pi/180,LAT(indo)*pi/180);
  [x,y]=deal(LON(indo),LAT(indo));
  p1=plot(x,y,'o');
  xl(1)=xlabel('longitude');
  yl(1)=ylabel('latitude');
  title('Original')

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  axes(ah(2))
  % Plot the newly interpolated function
  plotplm(rp,lon*pi/180,lat*pi/180,4)
  % Plot the sample locations
  hold on
  % [x,y]=mollweide(LON(indo)*pi/180,LAT(indo)*pi/180);
  [x,y]=deal(LON(indo),LAT(indo));
  p2=plot(x,y,'o');
  xl(2)=xlabel('longitude');
  title('Interpolated')

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  axes(ah(3))
  % Plot the selection of data
  plot(r(indo),'o-'); hold on
  % And their exact fit by interpolation
  plot(rp(indo),'r+-'); hold on
  axis tight
  title('Original and interpolated values [sampled]')

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  axes(ah(4))
  % Plot the entire data set
  plot(r(:),'-'); hold on
  % And their inexact approximation by interpolation
  plot(rp(:),'r-'); hold on
  axis tight
  title('Original and interpolated values [all]')
    
  % Cosmetics
  longticks(ah)
  set(ah(1:2),'clim',minmax([get(ah(1),'clim') get(ah(2),'clim')]))
  shrink(ah(3:4),1,1.25)
  nolabels(ha(3:4),2)
  set(ah(3),'xtick',unique([1 get(ah(3),'xtick')]))
  set(ah(4),'xtick',unique([1 get(ah(4),'xtick')]))
  set(ah(3:4),'ylim',ylix)
  set([p1 p2],'MarkerE','k','MarkerF','k','MarkerS',2)
  movev(ah(1:2),-.1)
  supertit(ah(1:2),...
	   sprintf('Abel-Poisson interpolation with h = %4.2f',h))
  fig2print(gcf,'landscape')
elseif strcmp(flatlon,'demo2')
  % We make up a function, sample randomly from it, and see how it
  % approximates the underlying function while fitting the sampled points
  % exactly 
  lmax=20;
  [m,l,mzero]=addmon(lmax);
  c=randn(addmup(lmax),2).*([l l].^(-1));
  c(1)=0; c(mzero,2)=0; lmcosi=[l m c];
  [r,lon,lat]=plm2xyz(lmcosi,1);
  tol=length(lon)*length(lat);
  defval('degres',0.4)
  fra=degres;
  unform=2;
  [LON,LAT]=meshgrid(lon,lat);
  % Really uniform on the sphere
  [lonr,latr]=randsphere(50);%ceil(fra*tol));
  indo=sub2ind(size(r),ceil(scale(latr,[1 length(lat)])),...
		 ceil(scale(lonr,[1 length(lon)])));
  latr=LAT(indo);
  
  % Decide on the width of an Abel-Poisson kernel by choosing h
  h=rand(3);
  h=[0.50 0.65 0.80];

  % And we judge visually how we've smoothed around the unknown locations
  clf
  [ah,ha]=krijetem(subnum(2,2));

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  axes(ah(1))
  % Plot the original function in Cartesian projection
  plotplm(r,lon*pi/180,lat*pi/180,4)
  % Plot the sample locations
  hold on
  % [x,y]=mollweide(LON(indo)*pi/180,LAT(indo)*pi/180);
  [x,y]=deal(LON(indo),LAT(indo));
  p(1)=plot(x,y,'o');
  xl(1)=xlabel('longitude');
  yl(1)=ylabel('latitude');
  title('Original')

  % Do the different kinds of interpolation
  for in=1:length(h)
    % Now use the function here for r(indo) and reconstruct it on a whole
    % fine grid... from lon(INDO) back to lon
    flatlonp=xyz2spl(r(indo),LAT(indo),LON(indo),LAT(:),LON(:),...
		     [],h(in));
    rp=reshape(flatlonp,length(lat),length(lon));
    
    % Now we verify that the ones that we provided are fit EXACTLY
    difer(flatlonp(indo)-r(indo))
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    axes(ah(1+in))
    % Plot the newly interpolated function
    plotplm(rp,lon*pi/180,lat*pi/180,4)
    % Plot the sample locations
    hold on
    % [x,y]=mollweide(LON(indo)*pi/180,LAT(indo)*pi/180);
    [x,y]=deal(LON(indo),LAT(indo));
    p(in+1)=plot(x,y,'o');
    xl(1+in)=xlabel('longitude');
    yl(1+in)=ylabel('latitude');
    title(sprintf('Abel-Poisson interpolation with h = %4.2f',h(in)))
  end

  % Cosmetics
  longticks(ah)
  % seemax(ah,3)
  set(ah,'clim',minmax(r(:)))
  nolabels(ha(3:4),2)
  delete(yl([2 4]))
  set(p,'MarkerE','k','MarkerF','w','MarkerS',4)
  movev(ah(1:2),-.125)
  movev(ah,.1)
  fig2print(gcf,'landscape')
  
  cb=colorbarf('hor',12,'Helvetica',[0.5703 0.1100 0.3347 0.0256]);
  layout(cb,0.5,'m','x')
  movev(cb,0.05); longticks(cb);
  kelicol
end


