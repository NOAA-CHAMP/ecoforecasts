function hh=featherdt_ERROR(varargin)
%FAILED EXPERIMENT - DO NOT RUN!
%FEATHERDT Feather plot with time axis.
%   FEATHER(DTS,U,V) plots the velocity vectors with components U and V as
%   arrows emanating from points along a horizontal time axis DTS.
%   FEATHERDT is useful for displaying direction and magnitude data that
%   is collected along a path.
%
%   FEATHERDT(DTS,Z) for complex Z is the same as FEATHER(DTS,REAL(Z),IMAG(Z)).
%   FEATHERDT(...,'LineSpec') uses the color and linestyle specification
%   from 'LineSpec' (see PLOT for possibilities).
%
%   FEATHERDT(AX,...) plots into AX instead of GCA.
%
%   H = FEATHERDT(...) returns a vector of line handles.
%
%   Example:
%      theta = (-90:10:90)*pi/180; r = 2*ones(size(theta));
%      dts = datenum(2008,1,1) + (1:length(theta));
%      [u,v] = pol2cart(theta,r);
%      featherdt(dts,u,v), axis equal
%
%   See also FEATHER, COMPASS, ROSE, QUIVER.

%   BASED ENTIRELY ON FEATHER:
%   Charles R. Denham, MathWorks 3-20-89
%   Modified 1-2-92, ls.
%   Modified 12-7-93 Mark W. Reichelt
%   Copyright 1984-2005 The MathWorks, Inc.
%   $Revision: 5.13.4.4 $  $Date: 2005/04/28 19:56:27 $

error('This is a failed experiment! Do not run...');

% Parse possible Axes input
[cax,args,nargs] = axescheck(varargin{:});
error(nargchk(1,4,nargs,'struct'));

if nargs > 0, dts = args{1}; end
if nargs > 1, x = args{2}; end
if nargs > 2, y = args{3}; end
if nargs > 3, s = args{4}; end

if ischar(x)
    error(id('SecondNumericInput'),'Second argument must be numeric.');
end
xx = [0 1 .8 1 .8]';
yy = [0 0 .08 0 -.08].';
arrow = xx + yy.*sqrt(-1);

if nargs == 3
   if ischar(y)
      s = y;
      y = imag(x); x = real(x);
   else
      s = '-';
   end
elseif nargs == 2
   s = '-';
   y = imag(x); x = real(x);
end
if ischar(x) || ischar(y)
    error(id('LeadingNumericInputs'),...
          'First 1 or 2 numeric arguments must be numeric.')
end
[st,co,mark,msg] = colstyle(s); error(msg) %#ok

x = x(:);
y = y(:);
if length(x) ~= length(y)
   error(id('LengthMismatch'),'X and Y must be the same length.');
end
m = size(x,1);

z = (x + y.*sqrt(-1)).';
%a = arrow * z + ones(5,1)*(1:m);
a = arrow * z + ones(5,1)*(dts');

% Plot the feather
cax = newplot(cax);
h = plot('v6', real(a), imag(a), [st co mark], ...
    'parent',cax);
if isempty(co),
  co = get(cax,'colororder');
  set(h,'color',co(1,:))
end

if nargout>0, hh = h; end

function str=id(str)
str = ['MATLAB:feather:' str];

