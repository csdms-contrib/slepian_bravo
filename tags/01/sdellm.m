function sdellm 
% SDELLM
%
% Calculates the error covariance matrix of a spherically symmetric
% (double-)polar gap satellite recovery problem.
%
% INPUT:
% 
% TH       Radius of the polar gap, in degrees
% L        Bandwidth
% N        White noise level
% S        White source level
% a        Satellite altitude, as a fraction of the radius
% m        Order
% 
% Last modified by fjsimons-at-alum.mit.edu, 08.08.2005

[E,V,th,C,ngl1,ngl2,K]=sdwcap2(TH,L,m,0);



