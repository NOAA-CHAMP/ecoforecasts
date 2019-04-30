function maximize_graph(fh)
%function maximize_graph(fh)
%
% Enlarge plot figure FH to fill available screen. If FH not specified,
% enlarge the GCF (see). Sets 'units','normalized', then 'position'.
% Position coordinates were found by actually maximizing a FIGURE on the
% screen of laptop MANANNAN.aoml.noaa.gov. 
% 
% Last Saved Time-stamp: <Mon 2018-09-17 23:01:04 Eastern Daylight Time gramer>

  if ( ~exist('fh','var') || isempty(fh) )
    fh = gcf;
  end;

  set(fh, 'units','normalized');
  %set(fh, 'position',[0.00 0.00 1.00 0.91]);
  set(fh, 'position',[0.00 0.00 0.99 0.91]);

return;
