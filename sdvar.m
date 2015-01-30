function V=sdvar(N,eta,Vo)
% V=SDVAR(N,eta,Vo)
%
% Calculates the VARIANCE MATRIX for geodetic inversions
%
% INPUT:
%
% N           Noise level
% eta         Damping parameter (0 <= eta <= 1) or NaN
%             If NaN this is the Slepian approach!
% Vo          (L+1)^2 vector with eigenvalues from GALPHA
%
% OUTPUT:
%
% V           (L+1)^2 X (L+1)^2 matrix with variance terms
%                 K   X    K    matrix in the Slepian case
%
% Last modified by fjsimons-at-alum.mit.edu, 17.10.2005

defval('N',0.1)
defval('eta',0.5)

if ~isnan(eta)
  % Damped spherical harmonic approach
  % Construct the help matrix
  H=Vo*(1-eta)+eta;
  
  % Construct the variance matrix
  V=N*diag(H.^(-2).*Vo);
else
  % Truncated Slepian function approach
  V=N*diag(1./Vo);
end



