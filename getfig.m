function fh = getfig(hdl)
%function fh = getfig(hdl)
%
% Return the handle of the figure that contains drawing object HDL (where HDL
% can be for example an AXES handle, a LINE object, TEXT object, COLORBAR, etc.
% 
% Last Saved Time-stamp: <Tue 2010-01-26 13:44:50 Eastern Standard Time Lew.Gramer>

  fh = [];

  if ( ishandle(hdl) )
    par = get(hdl, 'Parent');
    while ( ~isempty(par) && (par ~= 0) )
      hdl = par;
      par = get(hdl, 'Parent');
    end;
    if ( par == 0 )
      fh = hdl;
    end;
  end;

return;
