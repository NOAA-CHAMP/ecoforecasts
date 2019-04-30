function lpath = get_relative_path(ldir)
%function lpath = get_relative_path(ldir)
%
% Store/retrieve data in the local directory of OUR CALLER (or of PWD, v.)
%
% Last Saved Time-stamp: <Mon 2016-02-01 16:05:59 Eastern Standard Time gramer>

  %[pathroot, ig, ig, ig] = evalin('caller',['fileparts(mfilename(''fullpath''))']);
  % EVALIN did not perform as I expected here, so...
  [st,ix] = dbstack('-completenames');
  if ( ix > 1 )
    %[pathroot, ig, ig, ig] = fileparts(st(ix-1).file);
    [pathroot, ig, ig] = fileparts(st(ix-1).file);
  else
    % Caller is apparently the base workspace
    pathroot = pwd;
  end;
  clear st ix ig

  if ( ~exist('ldir','var') || isempty(ldir) )
    ldir = '';
  end;

  if ( length(ldir) > 1 && strcmp(ldir(1:2), '..') )
    % This implements the "../" in a portable way...
    %[pathroot, ig, ig, ig] = fileparts(pathroot);
    [pathroot, ig, ig] = fileparts(pathroot);
    if ( length(ldir) == 2 )
      ldir = '';
    else
      ldir = ldir(3:end);
    end;
  end;

  lpath = fullfile(pathroot, ldir, '');

return;
