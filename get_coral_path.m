function lpath = get_coral_path(ldir)
%function lpath = get_coral_path(ldir)
%
% Store/retrieve data in the local CORAL projects directory

  if ( ~exist('ldir','var') || isempty(ldir) )
    ldir = '';
  end;
  ldir = regexprep(ldir,'[\\/]',filesep);

  % % THIS DOESN'T WORK FOR ALL INSTALLS!
  % [pathroot, ig, ig] = fileparts(mfilename('fullpath'));
  % [pathroot, ig, ig] = fileparts(pathroot);
  % [pathroot, ig, ig] = fileparts(pathroot);
  % coralpath = fullfile(pathroot,'coral');
  % lpath = fullfile(coralpath,ldir);

  [coralpath, ig, ig] = fileparts(which('noaa_aoml_coral_project_path'));
  lpath = fullfile(coralpath,ldir);

return;
