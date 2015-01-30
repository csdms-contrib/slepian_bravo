function lmcosi=klm2lmcosi(C,wot)
% lmcosi=KLM2LMCOSI(C,wot)
%
% From the eigenfunctions C of the output of KERNELC, extracts the wot'th
% column and produces a matrix of the form lmcosi suitable for input to
% PLM2XYZ or PLOTPLM etc. Since KERNELC works with unit-normalized
% harmonics and PLM2XYZ with 4*pi-normalized ones, we renormalize. We
% didn't use to do this elsewhere, but from now on, it will matter.
%
% SEE ALSO: PLOTSLEP, LOCALIZATION, GLM2LMCOSI
% 
% Last modified by fjsimons-at-alum.mit.edu, 07/11/2012

% Collect the output into a format that PLM2XYZ knows how to interpret
[dems,~,~,lmcosi,mzin]=addmon(sqrt(size(C,1))-1);
% Construct the full matrix
lmcosi(:,3:4)=reshape(insert(C(:,wot),0,mzin),2,length(dems))';
% Renormalize so that when multiplied with 4pi-harmonics they have the
% same integral. Note that this has nothing special for the m=0 zonal
% coefficients - it's simply the last mile of the normalization. Check
% this using Fibonacci_grid as suggested in GALPHA.
lmcosi(:,3:4)=lmcosi(:,3:4)/sqrt(4*pi);

% Remember that if you were to use YLM for rendering you'd have to add
% the (-1)^m phase factor which is missing from this. 
