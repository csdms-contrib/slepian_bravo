function [err,F,G,N,NA,th]=sdsumk(ZTH,ZL,nth,cb,xver)
% [err,F,G,N,NA,th]=SDSUMK(TH,L,nth,cb,xver)
%
% Sums eigenvalue-weighted squared eigenfunctions on the
% DOUBLE-POLAR CAP or the complementary LATITUDINAL BELT.
% Sums Shannon number terms as well as all of the terms.
% Also compares the full unweighted sum to N/A.
%
% INPUT:
%
% TH          Cap radius, in degrees (scalar)
% L           Bandwidth (scalar)
% nth         Number of spatial sampling points
% cb          1 Double-cap, cap eigenvalue weighted
%             2 Latitudinal belt, belt eigenvalue weighted
%             3 Double-cap, (1-cap eigenvalue) weighted
%             4 Latitudinal belt, (1-belt eigenvalue) weighted
% xver        0 No excessive validation & verification (default)
%             1 With excessive validation & verification
%
% OUTPUT:
%
% F           The eigenvalue-weighted sum to the Shannon number
% G           The eigenvalue-weighted sum over all functions
% N           The Shannon number that ends up being used
% th          The colatitudes at which F and G are evaluated
%
% Finally, a bit of a professional routine.
% Should replace SDSUMALL, SDSUMALL2
%
% Last modified by fjsimons-at-alum.mit.edu, 16.08.2005

defval('ZTH',3); 
defval('ZL', 45);
defval('nth',720);
defval('cb',1);
defval('xver',0);

% Define tolerance
tol=1e-10;

for ondex=1:length(ZTH)
  TH=ZTH(ondex);
  L=ZL(ondex);
  % Calculate the Shannon number for the DOUBLE CAP
  Nor=(L+1)^2*(1-cos(TH*pi/180));
  switch cb
   case {2 4}
    % Calculate the Shannon number for the LATITUDINAL BELT
    Nor=(L+1)^2-Nor;
   case {1 3}
   otherwise
    error('Specify valid case!')
  end
  % Whatever it is, round it to the nearest integer
  N=round(Nor);
    
  % Calculate N/A, always the same for identical L
  NA=(L+1)^2/4/pi;
  
  disp(sprintf('L = %3.3i ; TH = %2.2i ; N = %3.3i ; N/A = %8.3f',...
	       L,TH,N,NA))

  % If we've done it before, don't redo it
  fnpl=sprintf('%s/SDSUMK-%i-%i-%i-%i-%i.mat',...
	       fullfile(getenv('IFILES'),'WIECZOREK'),...
	       TH,L,N,nth,cb);
    
  if exist(fnpl,'file')~=2 | 1==1
    % Figure out order to sum them in, get expanded eigenvalues
    % Both 1 and 3 map to 1, the caps, and 2 and 4 to 2, the belt  
    [lrnk,mrnk,lval,VV,Vsum]=sdelm(TH,L,mod(cb+1,2)+1,1);

    % If the sum is wrong, and we mean, really wrong, not due to rounding
    if round(Vsum)~=N & abs(Nor-Vsum)>tol
      error('Shannon number must equal sum of eigenvalues')
    end

    % Must have all the pairs of orders to be longitudinally independent
    if N~=0
      if mrnk(N)~=0 & ...  % Not zero
		  lrnk(N)==lrnk(N+1) & ... % Same degree
		  [mrnk(N)==-mrnk(N+1)] % Just switched signs
	N=N+1;
	disp(sprintf(...
	    'Shannon number updated to %i to avoid degenerate split',N))
      end
    end
    % Calculate all the caps tapers for all orders, at one longitude only
    % This is always the same and in the same order
    G=0; G2=0;
    for index=0:L
      [E{index+1},Vg,th,C,T,V{index+1}]=...
	  grunbaum2(TH,L,index,nth,1);
      % Calculate full unweighted sum
      G=G+sum(E{index+1}.^2,2);
      
      % Excessive verification starts here
      if xver==1
	% Reculculate using the standard way
	[E2{index+1},V2{index+1},th,C2,n1,n2,T2]=...
	    sdwcap2(TH,L,index,nth,-1,1);
	% Calculate full unweighted sum
	G2=G2+sum(E2{index+1}.^2,2);
	% The localization matrices should commute, do they?
	DT=mean(mean(abs(T*T2-T2*T)));
	if DT<tol
	  disp(sprintf('Commutation works to within   %8.3e',DT))
	else 
	  warning(sprintf('Commutation troubled: %8.3e > %1.0e',DT,tol))
	end

	% The eigenvalues should be nearly identical
	DV=mean(abs(V{index+1}-V2{index+1}));
	if DV<tol
	  disp(sprintf('Eigenvalues identical within  %8.3e',DV))
	else 
	  warning(sprintf('Eigenvalue troubled: %8.3e > %1.0e',DV,tol))
	end
	
	% The spectral eigenfunction matrices should be identical as long
        % as the eigenvalue is not too near zero; error if N1 is not N2
	N1=sum(V{index+1}>tol);
	N2=sum(V2{index+1}>tol);
	if N1>0
	  D=E{index+1}(:,1:N1)-E2{index+1}(:,1:N2);
	  DE=mean(mean(abs(D)));
	  if DE<tol
	    disp(sprintf(...
		'Significant eigenfunctions to %8.3e',DE))
	  else 
	    warning(sprintf(...
		'Significant eigenfunctions troubled: %8.3e > %1.0e',DE,tol))
	  end
	end
      end
      
      % Switch eigenvalue from cap to belt, 1-cap or 1-belt
      switch cb
       case {1 4}
	% Cap ordering/Cap eigenvalue &  Belt ordering/1-Belt eigenvalue
	V{index+1}=V{index+1};
       case {2 3}
	% Belt ordering/Belt eigenvalue & Cap ordering/1-Cap eigenvalue
	V{index+1}=1-V{index+1};
      end
    end % Sum over orders

    if xver==1
      % Make sure that the full unweighted sum reduces to N/A
      if any(abs(G2-NA)>tol)
	warning(sprintf(...
	    'Unweighted sum off N/A by %8.3f for L= %i and TH= %i',...
	    max(abs(G2-NA)),L,TH))
      else
	disp(sprintf(...
	    'Unweighted sum within N/A by %8.3e for L= %i and TH= %i',...
	    max(abs(G2-NA)),L,TH))
      end
    end
    % Make sure that the full unweighted sum reduces to N/A
    if any(abs(G-NA)>tol)
      warning(sprintf(...
	  'Unweighted sum off N/A by %8.3f for L= %i and TH= %i',...
	  max(abs(G-NA)),L,TH))
      err=mean(abs(G(:)-NA));
      return
    else
      disp(sprintf(...
	  'Unweighted sum within N/A by %8.3e for L= %i and TH= %i',...
	  max(abs(G-NA)),L,TH))
    end

    err=mean(abs(G(:)-NA));

    if N~=0
      % No point in keeping G here, is there, thus reassign
      % So we don't need to sum with the longitudinal part
      smrnk=mrnk(1:N); rmrnk=mrnk(N+1:end);
      slrnk=lrnk(1:N); rlrnk=lrnk(N+1:end);
      
      % And sum only the unique lrnk abs(m) combinations
      % But now these are sorted upwards
      sift=unique([slrnk abs(smrnk)],'rows');
      seft=unique([rlrnk abs(rmrnk)],'rows');
      slrnk=sift(:,1);    rlrnk=seft(:,1);
      smrnk=sift(:,2);    rmrnk=seft(:,2);
      
      % Calculate cumulative sums weighted by eigenvalue!
      F=0; GG=0;
      % Now calculate the sum only to the Shannon number N
      for index=1:length(smrnk)
	% Sum the eigenvalue-weighted square of the lrnk'th sequential
	% eigenfunctions of the mrnk'th order
	rnkm=smrnk(index)+1;
	rnkl=slrnk(index);
	
	% Sum cumulatively with the eigenvalue fixed above 
	F=F+(E{rnkm}(:,rnkl)).^2*V{rnkm}(rnkl);
	GG=GG+sum(E{rnkm}(:,rnkl).^2,2);
      end
      G=F;
      % And keep going all the way by adding the remaining ones
      % Really could have done this in the first loop, but who cares?
      for index=1:size(seft,1)
	rnkm=seft(index,2)+1;
	rnkl=seft(index,1);
	
	% Sum cumulatively with the eigenvalue fixed above 
	G=G+(E{rnkm}(:,rnkl)).^2*V{rnkm}(rnkl);
	GG=GG+sum(E{rnkm}(:,rnkl).^2,2);
      end
      eval(sprintf('save %s F G N th',fnpl))
    end
  else
    eval(sprintf('load %s',fnpl))
  end
end


