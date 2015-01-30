function [K,Al,Dlta,Kp]=abelpoisson(h,Dlta,L)
% [K,Al,Dlta,Kp]=ABELPOISSON(h,Dlta,L)
%
% Calculates the Abel-Poisson kernel and verifies the funky equation
% producing it --- see the work of Freeden and his posse
%
% INPUT:
% 
% h      The concentration parameter, can be a vector
% L      The bandwidth of the explicit expression - if desired
%
% OUTPUT:
%
% K      The kernels in terms of spherical distance, [length(Delta)*length(h)]
% Al     The degree-dependent sequence
% Dlta   The geodesic angular distance
% Kp     The kernel in a bandlimited (approximate) expansion
%
% EXAMPLE:
%
% [K,Al,Dlta,Kp]=abelpoisson;
%
% Last modified by fjsimons-at-alum.mit.edu, 08/06/2007

defval('h',0.50)
defval('L',50)

crsfin=unique([linspace(-pi/6,pi/6,100) linspace(-pi,pi,900)]);

defval('Dlta',crsfin)

% Make dimensions work out
N=length(Dlta);
Dlta=Dlta(:);
h=h(:)';
hh=repmat(h,N,1);

% Calculate the kernel
K=(1-hh.^2)./4/pi./(1+hh.^2-2*cos(Dlta)*h).^(3/2);

% Return the spectrum
Al=repmat(h,L+1,1).^(repmat(-[0:L]'/2,1,length(h)));

if nargout==4
  % Calculate an alternative representation of the kernel
  Kp=NaN(size(K));
  for index=1:length(h)
    Kp(:,index)=sum(repmat((2*[0:L]'+1)/4/pi.*Al(:,index).^-2,1,N).*...
		    plm(0:L,0,cos(Dlta)),1);
    % Report on the fit
    disp(sprintf('Relative misfit L = %i approximation at h = %3.2f : %5.3e %s',...
	    L,h(index),norm(K(:,index)-Kp(:,index))./norm(K(:,index))*100,'%'))
  end
end




