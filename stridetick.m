function stridetick(varargin)
%function stridetick([AX,AXIS,STRIDE]])
%
% For an AXES, STRIDETICK modifies the ticks of an axis to use every Nth
% value. E.g., STRIDETICK by itself defaults to GCA, AXIS 'x', STRIDE 2.
% STRIDETICK('y',3) - Remove every 2nd and 3rd value in YTick property.
% STRIDETICK(H,'z',4) - Remove every 2nd-4th value in ZTick of axes H.
% If STRIDE>=NUMEL(ticks), remove all axis ticks except first and last.
%
% NOTE: Effects of repeated calls to STRIDETICK are cumulative.
%
% Function may be useful when printing figures with crowded tick labels.
%
% Last Saved Time-stamp: <Sat 2015-07-25 09:25:14 Eastern Daylight Time gramer>

  if (nargin==3)
    ax = varargin{end-2};
    xyz = varargin{end-1};
    strd = varargin{end};
  elseif (nargin==2)
    ax = gca;
    xyz = varargin{end-1};
    strd = varargin{end};
  elseif (nargin==1)
    ax = gca;
    xyz = 'x';
    strd = varargin{end};
  elseif (nargin==0)
    ax = gca;
    xyz = 'x';
    strd = 2;
  else
    error('Too many input arguments');
  end;

  propname = [xyz,'Tick'];

  tick = get(ax,propname); 
  if ( numel(tick) <= strd )
    tick = tick([1,end]);
  else
    tick = tick(1:strd:end);
  end;
  set(ax,propname,tick);

return;
