function p = get_datapath(ldir)
%function p = get_datapath(ldir)
% Return the path to the data directory LDIR. DEFAULT: "~/Documents/data" on
% a Mac, "d:/data" on Lew's PC. EXAMPLE: Calling GET_DATAPATH('hycom/gom') or
% GET_DATAPATH('hycom\gom') on the PC would return 'D:/data/hycom/gom'. Call
% GET_DATAPATH('../archive/NWPS') would return the path to the NWPS archive.

  if ( ismac )
    pathroot = '~/Documents/data';
  else
    pathroot = 'd:/data';
  end;

  if ( ~exist('ldir','var') || isempty(ldir) )
    ldir = '';
  end;
  if ( length(ldir) > 1 && strcmp(ldir(1:2), '..') )
    % This implements "../" in a portable way...
    [pathroot, ig, ig] = fileparts(pathroot);
    ldir = ldir(3:end);
  end;

  p = fullfile(pathroot, ldir, '');
return;
