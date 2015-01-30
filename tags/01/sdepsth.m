function [e2th,theta]=sdepsth(TH,L,nth)
% [e2th,theta]=SDEPSTH(TH,L,nth)
%
% INPUT:
%
% TH          Angular extent of the spherical cap, in degrees
% L           Bandwidth (maximum angular degree)
% 
% OUTPUT:
%
% e2th        Mean-square error in space
% theta       Colatitudes where this applies
%
% Computes the mean-square spatial error of a global field recovered from
% observations suffering from a polar gap with a double polar cap.
% Watch out for problems with cos in calculation of Legendre.
%
% Note this is, Y^T*(S*Dbar)*Y but we make the infinity mistake;
% the difference is the power on the eigenvalue
%
% See also SDGELDEREN and SDVAR, SDBIAS, SDERR, SDERK
%
% Last modified by fjsimons-at-alum.mit.edu, 05/21/2009

defval('TH',30)
defval('L',18)
defval('nth',720)

% Use precomputed output 
fnpl0=sprintf('%s/SDEPSTH-%i-%i-%i.mat',...
	     fullfile(getenv('IFILES'),'SDWCAP2'),TH,L,nth);

if exist(fnpl0,'file')==2
  disp(sprintf('Loading %s',fnpl0))
  load(fnpl0)
else
  % Define data grid
  theta=linspace(0,pi,nth);
  x=cos(theta);
  as=1; % Equally spaced
  
  % Get Legendre function values if equally spaced on 0->pi
  fnpl=sprintf('%s/LSSM-%i-%i.mat',...
	       fullfile(getenv('IFILES'),'LEGENDRE'),L,nth);

  if exist(fnpl,'file')==2 && as==1 
    load(fnpl)
  else  
    libb=0;
    % Evaluate Legendre polynomials at selected points
    Plm=repmat(NaN,length(theta),addmup(L));
    if L>200
      h=waitbar(0,'Evaluating all Legendre polynomials');
    end
    in1=0;
    in2=1;
    % Always start from the beginning
    % BUT IT ALSO INCLUDES sqrt(2-dom) !
    for l=0:L
      if libb==0
	Plm(:,in1+1:in2)=(legendre(l,x(:)','sch')*sqrt(2*l+1))';
      else
	Plm(:,in1+1:in2)=(libbrecht(l,x(:)','sch')*sqrt(2*l+1))';
      end
      in1=in2;
      in2=in1+l+2;
      if L>200
	waitbar((l+1)/(L+1),h)
      end
    end
    if L>200
      delete(h)
    end
    if as==1
      save(fnpl,'Plm')
    end
  end

  % Make these equal to the Xlm in our notation
  % But note the factor of sqrt(2) is already in here
  Plm=Plm/sqrt(4*pi);

  % Initialize mean-square error array
  e2th=repmat(0,nth,1);

  % Figure out where the orders start and end
  [ms,ls,mz]=addmon(L);

  % Loop over the orders and sum over the degrees
  for m=0:L
    % Get localization matrix K at fixed-order
    % Ask not a single theta of the SDWCAP2 function
    [E,V,th,C,ngl1,ngl2,K]=sdwcap2(TH,L,m,0);
    % Pick out the polynomials at the right order
    Pm=Plm(:,mz(m+1:end)+m);
    % This should be symmetric, people
    % Note that the factor of two is already in there
    e2th=e2th+diag(Pm*K*Pm');
  end
  save(fnpl0,'e2th','theta')
end



