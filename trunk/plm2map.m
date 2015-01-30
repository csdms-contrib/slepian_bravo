function imap=plm2map(lmcosi)
% imap=PLM2MAP(lmcosi)
%
% INPUT:
%
% lmcosi    The usual coefficient matrix
%
% OUTPUT:
%
% imap      A matrix suitable for plotting
%
% Puts a matrix with real spherical harmonic expansion coefficients into
% a matrix for easy "eyeballing".
%
% SEE ALSO:
%
% SDSNEEUW, SDEPSLM, SLEP2MAP
%
% Last modified by fjsimons-at-alum.mit.edu, 05/15/2009
% Last modified by kwlewis-at-princeton.mit.edu, 12/08/2010

lmin=lmcosi(1,1);
lmax=lmcosi(end,1);

imap=nan(lmax+1,2*lmax+1);

missl=addmup(lmcosi(1,1)-1);

for degr=lmin:lmax
  reldeg=[addmup(degr-1)+1:addmup(degr)]-missl;
  imap(degr+1,1+lmax:-1:1+lmax-degr)=lmcosi(reldeg,3)';  
  imap(degr+1,1+lmax+1:1+lmax+degr)=lmcosi(reldeg(2:end),4)';
end


