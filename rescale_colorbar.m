function cbh = rescale_colorbar(cbh,sclfun,fmt)
%function cbh = rescale_colorbar(cbh,sclfun,fmt)
%
% Modify displayed label strings in colorbar with handle CBH (DEFAULT: CBH =
% COLORBAR) using scaling function SCLFUN (DEFAULT: POWER(10,x)) and NUM2STR
% format specifier FMT (DEFAULT: '%5.2f').
%
% NOTE: Default is a colorbar with logarithmic scale and 2-digit precision.
%
% Last Saved Time-stamp: <Sun 2017-04-23 21:35:22 Eastern Daylight Time gramer>

  if ( ~exist('cbh','var') || isempty(cbh) )
    cbh = get(gca,'Colorbar');
  end;
  if ( ~exist('sclfun','var') || isempty(sclfun) )
    sclfun = @(x)(power(10,x));
  end;
  if ( ~exist('fmt','var') || isempty(fmt) )
    fmt = '%5.2f';
  end;

  th = get(cbh,'TickLabels');
  for thix=1:numel(th);
    th{thix} = num2str(sclfun(str2num(th{thix})),fmt);
  end;
  set(cbh,'TickLabels',th);

return;
