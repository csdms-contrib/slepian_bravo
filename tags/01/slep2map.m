function imap=slep2map(falpha,MTAP)
% imap=SLEP2MAP(falpha,MTAP)
%
% INPUT:
%
% falpha    The usual coefficient matrix
% MTAP      The block-sorted orders that this belongs to
%
% OUTPUT:
%
% imap      A matrix suitable for plotting
%
% Puts a matrix with Slepian expansion coefficients into
% a matrix for easy "eyeballing".
%
% SEE ALSO:
%
% SDSNEEUW, SDEPSLM, PLM2MAP
%
% Last modified by fjsimons-at-alum.mit.edu, 05/19/2009

% Figure out bandwidth and do basic checking
L=sqrt(length(falpha))-1;
difer(minmax(MTAP)-[-L L],[],[],NaN)
difer([1:L]-unique(abs(MTAP(~~MTAP))),[],[],NaN)

imap=nan(L+1,2*L+1);

% This figures out the ranges into falpha for a particular positive or
% negative order
relems=[L+1-sort(abs(-L:L))];
relems=[cumsum([1 relems(1:end-1)])' cumsum(relems)'];

for m=0:L
  if ~~m
    % The negative orders
    relemm=relems(2*m,1):relems(2*m,2);
    difer(unique(MTAP(relemm))+m,[],[],NaN)
    % Stick the coefficients in at the right place
    imap(m+1:L+1,L+1-m)=falpha(relemm);
  end
  % The positive orders
  relemp=relems(2*m+1,1):relems(2*m+1,2);
  difer(unique(MTAP(relemp))-m,[],[],NaN)
  % Stick the coefficients in at the right place
  imap(m+1:L+1,L+1+m)=falpha(relemp);
end

