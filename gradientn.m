function varargout = gradientn(f,npts,varargin)
%GRADIENTN: Lew.Gramer@noaa.gov 2011 Mar 24: Identical to GRADIENT (R2007a
%  version), except allows centered finite differences over NPTS points, not
%  just three. NPTS may be any of 3 (behavior similar to GRADIENT),5,7,9,11.
% Lew.Gramer@noaa.gov, 2016 Jul 14:
%  The formulae used for 5, 7, and 9-point templates were based on:
%   http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/central-differences/
% Lew.Gramer@noaa.gov, 2018 Jun 23:
%  The formula used for 11-point template was based on:
%   https://samson.chem.umass.edu/pub_pdf/pap73.pdf
%
%GRADIENT Approximate gradient.
%       
%   [FX,FY] = GRADIENTN(F,NPTS) returns NPTS-point numerical gradient of the
%   matrix F. FX corresponds to dF/dx, the differences in x (horizontal) 
%   direction. FY corresponds to dF/dy, the differences in y (vertical) 
%   direction. The spacing between points in each direction is assumed to 
%   be one. When F is a vector, DF=GRADIENTN(F,NPTS) is the 1-D gradient.
%
%   [FX,FY] = GRADIENTN(F,NPTS,H), where H is a scalar, uses H as the
%   spacing between points in each direction.
%
%   [FX,FY] = GRADIENTN(F,NPTS,HX,HY), when F is 2-D, uses the spacing
%   specified by HX and HY. HX and HY can either be scalars to specify
%   the spacing between coordinates or vectors to specify the
%   coordinates of the points.  If HX and HY are vectors, their length
%   must match the corresponding dimension of F.
%
%   [FX,FY,FZ] = GRADIENTN(F,NPTS), when F is a 3-D array, returns the
%   numerical gradient of F. FZ corresponds to dF/dz, the differences
%   in the z direction. GRADIENTN(F,NPTS,H), where H is a scalar, 
%   uses H as the spacing between points in each direction.
%
%   [FX,FY,FZ] = GRADIENTN(F,NPTS,HX,HY,HZ) uses the spacing given by
%   HX, HY, HZ. 
%
%   [FX,FY,FZ,...] = GRADIENTN(F,NPTS,...) extends similarly when F is N-D
%   and must be invoked with N outputs and either 2 or N+1 inputs.
%
%   Note: The first output FX is always the gradient along the 2nd
%   dimension of F, going across columns.  The second output FY is always
%   the gradient along the 1st dimension of F, going across rows.  For the
%   third output FZ and the outputs that follow, the Nth output is the
%   gradient along the Nth dimension of F.
%
%   Examples:
%       [x,y] = meshgrid(-2:.2:2, -2:.2:2);
%       z = x .* exp(-x.^2 - y.^2);
%       [px,py] = gradientn(z,3,.2,.2);
%       contour(z), hold on, quiver(px,py), hold off
%
%   Class support for input F:
%      float: double, single
%
%   See also DIFF, DEL2.
%
% Last Saved Time-stamp: <Sat 2018-06-23 17:33:59 Eastern Daylight Time gramer>

%   ORIGINAL GRADIENT CODE: Copyright 1984-2006 The MathWorks, Inc.
%   ORIGINAL GRADIENT CODE: $Revision: 5.17.4.5 $  $Date: 2006/06/20 20:09:24 $

if ( nargin<2 )
  error('Ecoforecasts:gradientn:InvalidInputs','GRADIENTN requires at least two arguments!');
end;
if ( ~isnumeric(npts) || ~isscalar(npts) )
  error('Ecoforecasts:gradientn:InvalidInputs','GRADIENTN second arg must be a scalar whole number!');
end;
if ( ~ismember(npts,[3,5,7,9,11]) )
  error('Ecoforecasts:gradientn:InvalidInputs','Do not know how to finite difference over %d points!', npts);
end;  

[msg,f,ndim,loc,rflag] = parse_inputs(f,varargin);
if ~isempty(msg), error('Ecoforecasts:gradientn:InvalidInputs', msg); end

% Loop over each dimension. Permute so that the gradient is always taken along
% the columns.

if ndim == 1
  perm = [1 2];
else
  perm = [2:ndim 1]; % Cyclic permutation
end

for k = 1:ndim
   [n,p] = size(f);
   h = loc{k}(:);   
   g  = zeros(size(f),class(f)); % case of singleton dimension

   % Take forward differences on edges
   if n > 1
      g(1,:) = (f(2,:) - f(1,:))/(h(2)-h(1));
      g(n,:) = (f(n,:) - f(n-1,:))/(h(n)-h(n-1));
   end

   % Take reduced-template differences when too close to edges
   if n > 2
      switch ( npts )
       case {5,7,9,11},
        g(2,:) = (f(3,:) - f(1,:))/(h(3)-h(1));
        g(n-1,:) = (f(n,:) - f(n-2,:))/(h(n)-h(n-2));
      end;
   end
   if n > 3
     switch ( npts )
      case {7,9,11},
       g(3,:) = (f(1,:) - 8*f(2,:) + 8*f(4,:) - f(5,:)) ./ ((h(5) - h(1))*3);
       g(n-2,:) = (f(n-4,:) - 8*f(n-3,:) + 8*f(n-1,:) - f(n,:)) ./ ((h(n) - h(n-4))*3);
     end;
   end
   if n > 4
     switch ( npts )
      case {9,11},
       g(4,:) = (-f(1,:) + 9*f(2,:) - 45*f(3,:) + 45*f(5,:) - 9*f(6,:) + f(7,:)) ./ ((h(7) - h(1))*10);
       g(n-3,:) = (-f(n-6,:) + 9*f(n-5,:) - 45*f(n-4,:) + 45*f(n-2,:) - 9*f(n-1,:) + f(n,:)) ./ ((h(n) - h(n-6))*10);
     end;
   end
   if n > 5
     switch ( npts )
      case {11},
       g(5,:) = (3*f(1,:) - 32*f(2,:) + 168*f(3,:) - 672*f(5,:) + 672*f(6,:) - 168*f(7,:) + 32*f(8,:) - 3*f(9,:)) ./ ((h(9) - h(1))*105);
       g(n-4,:) = (3*f(n-8,:) - 32*f(n-7,:) + 168*f(n-6,:) - 672*f(n-5,:) + 672*f(n-4,:) - 168*f(n-2,:) + 32*f(n-1,:) - 3*f(n,:)) ./ ((h(n) - h(n-8))*105);
     end;
   end

   % Take centered differences on interior points
   if n > ceil(npts/2)
     switch ( npts )
      case 3,
       h = h(3:n) - h(1:n-2); % del_h*2
       g(2:n-1,:) = (f(3:n,:)-f(1:n-2,:))./h(:,ones(p,1));
      case 5,
       h = (h(5:n) - h(1:n-4))*3; % del_h*12
       g(3:n-2,:) = (f(1:n-4,:) - 8*f(2:n-3,:) + 8*f(4:n-1,:) - f(5:n,:)) ./ h(:,ones(p,1));
      case 7,
       h = (h(7:n) - h(1:n-6))*10; % del_h*60
       g(4:n-3,:) = (-f(1:n-6,:) + 9*f(2:n-5,:) - 45*f(3:n-4,:) + 45*f(5:n-2,:) - 9*f(6:n-1,:) + f(7:n,:)) ./ h(:,ones(p,1));
      case 9,
       h = (h(9:n) - h(1:n-8))*105; % del_h*840
       g(5:n-4,:) = (3*f(1:n-8,:) - 32*f(2:n-7,:) + 168*f(3:n-6,:) - 672*f(4:n-5,:) + 672*f(6:n-3,:) - 168*f(7:n-2,:) + 32*f(8:n-1,:) - 3*f(9:n,:)) ./ h(:,ones(p,1));

      case 11,
       h = (h(11:n) - h(1:n-10))*126; % del_h*1260
       g(6:n-5,:) = (-1*f(1:n-10,:) + 12.5*f(2:n-9,:) - 75*f(3:n-8,:) + 300*f(4:n-7,:) - 1050*f(5:n-6,:) + 1050*f(7:n-4,:) - 300*f(8:n-3,:) + 75*f(9:n-2,:) - 12.4*f(10:n-1,:) + 1*f(11:n,:)) ./ h(:,ones(p,1));
     end;
   end

   varargout{k} = ipermute(g,[k:ndims(f) 1:k-1]);

   % Set up for next pass through the loop
   f = permute(f,perm);
end 

% Swap 1 and 2 since x is the second dimension and y is the first.
if ndim>1
  tmp = varargout{1};
  varargout{1} = varargout{2};
  varargout{2} = tmp;
end

if rflag, varargout{1} = varargout{1}.'; end


%-------------------------------------------------------
function [msg,f,ndim,loc,rflag] = parse_inputs(f,v)
%PARSE_INPUTS
%   [MSG,F,LOC,RFLAG] = PARSE_INPUTS(F,V) returns the spacing
%   LOC along the x,y,z,... directions and a row vector
%   flag RFLAG. MSG will be non-empty if there is an error.

msg = '';
loc = {};
nin = length(v)+1;

% Flag vector case and row vector case.
ndim = ndims(f);
vflag = 0; rflag = 0;
if ndims(f) == 2
   if size(f,2) == 1
      ndim = 1; vflag = 1; 
   elseif size(f,1) == 1    % Treat row vector as a column vector
      ndim = 1; vflag = 1; rflag = 1;
      f = f.';
   end;
end;
   
indx = size(f);

% Default step sizes: hx = hy = hz = 1
if nin == 1, % gradientn(f)
   for k = 1:ndims(f)
      loc(k) = {1:indx(k)};
   end;

elseif (nin == 2) % gradientn(f,h)
   % Expand scalar step size
   if (length(v{1})==1)
      for k = 1:ndims(f)
         h = v{1};
         loc(k) = {h*(1:indx(k))};
      end;
   % Check for vector case
   elseif vflag
      loc(1) = v(1);
   else
      msg = 'Invalid inputs to GRADIENTN.';
   end

elseif ndims(f) == numel(v), % gradientn(f,hx,hy,hz,...)
   % Swap 1 and 2 since x is the second dimension and y is the first.
   loc = v;
   if ndim>1
     tmp = loc{1};
     loc{1} = loc{2};
     loc{2} = tmp;
   end

   % replace any scalar step-size with corresponding position vector
   for k = 1:ndims(f)
      if length(loc{k})==1
         loc{k} = loc{k}*(1:indx(k));
      end;
   end;

else
   msg = 'Invalid inputs to GRADIENTN.';

end
