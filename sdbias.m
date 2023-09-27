function B=sdbias(S,eta,Vo)
% B=sdbias(S,eta,Vo)
%
% Calculates the BIAS MATRIX for geodetic spherical-harmonics and
% Slepian-function inversions following
% Simons and Dahlen (2006)
%
% INPUT:
%
% S           Signal level ("bandlimited" whitish)
% eta         Damping parameter (0 <= eta <= 1) or NaN
%             If NaN this is the Slepian approach!
% Vo          (L+1)^2 vector with eigenvalues from GALPHA
%
% OUTPUT:
%
% B           (L+1)^2 X (L+1)^2 matrix with bias terms
%                 K   X    K    matrix in the Slepian case
%
% Last modified by fjsimons-at-alum.mit.edu, 09/27/2023

defval('S',1)
defval('eta',1)

if ~isnan(eta)
  % Damped spherical harmonic approach
  % Construct the help matrix
  H=Vo*(1-eta)+eta;

  % Construct the bias matrix
  B=S*diag(H.^(-2).*(Vo-1).^2);
else
  % Truncated Slepian function approach
  B=S*eye(size(diag(Vo)));
end
