function result = read_hobo_csv(fnamepatt,readPressure,readLight)
%function result = read_hobo_csv(fnamepatt,readPressure,readLight)
%
% Read CSV data file for Hobo sea temperature (or pressure+temperature or
% temperature+light) recorder matching path specifier FNAMEPATT (passed to
% DIR, e.g., 'moorings/Data/T?11*csv'). If READPRESSURE is present and True,
% or if FNAMEPATT contains string 'pressure' (aNy CaSe), read pressure as
% first data column, sea temperature second. If READLIGHT present and True,
% or if FNAMEPATT contains string 'light', read temp as first column, then
% light. Converts pressure if present from PSI to [dbar], and light if
% present from lum/ft^2 to (approx.) PAR units [micro-mol quanta/m^2/s].
%
% Last Saved Time-stamp: <Fri 2011-09-09 10:13:06  lew.gramer>

  if ( ~exist('readPressure','var') || isempty(readPressure) )
    readPressure = false;
  end;
  if ( ~exist('readLight','var') || isempty(readLight) )
    readLight = false;
  end;

  result.seatemp = [];
  result.seapres = [];
  result.par = [];

  fs = dir(fnamepatt);
  if ( numel(fs) < 1 )
    error('No file matching pattern "%s"',fnamepatt);
  elseif ( numel(fs) > 1 )
    error('Multiple matches for pattern "%s"',fnamepatt);
  end;

  [fpath,ig] = fileparts(fnamepatt);
  fname = fullfile(fpath,fs.name);

  fid = fopen(fname,'r');
  if ( fid < 2 )
    error('Cannot open for reading "%s"',fname);
  end;
  C = textscan(fid,'%d,%[^,],%s','Headerlines',2,'Delimiter','');
  fclose(fid);

  if ( ~isempty(regexpi(fs.name,'pressure')) || readPressure )
    result.seatemp.date = datenum(C{2});
    result.seapres.date = result.seatemp.date;
    [result.seapres.data,result.seatemp.data] = cellfun(@read_hobo_csv_cvt_pres_temp,C{3});

    % Convert psi -> dbar
    result.seapres.data = result.seapres.data .* 0.68947573;

  elseif ( ~isempty(regexpi(fs.name,'light')) || readLight )
    result.seatemp.date = datenum(C{2});
    result.par.date = result.seatemp.date;
    [result.seatemp.data,result.par.data] = cellfun(@read_hobo_csv_cvt_temp_par,C{3});

    % Convert lumens/ft^2 -> umol/m^2/s
    % Uncertain??? Onset site suggests this rule-of-thumb conversion
    result.par.data = result.par.data ./ 13.03;

  else
    result.seatemp.date = datenum(C{2});
    result.seatemp.data = cellfun(@read_hobo_csv_cvt_temp,C{3});

  end;

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRIVATE FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Idiotic TEXTSCAN and idiotic STR2NUM!

function t = read_hobo_csv_cvt_temp(str)
  t = sscanf(str,'%f');
  if ( isempty(t) ); t = NaN; end;
return;

function [t,l] = read_hobo_csv_cvt_temp_par(str)
  vals = sscanf(str,'%f,%f');
  t = NaN;
  l = NaN;
  if ( numel(vals) > 0 )
    t = vals(1);
  end;
  if ( numel(vals) > 1 )
    l = vals(2);
  end;
return;

function [p,t] = read_hobo_csv_cvt_pres_temp(str)
  vals = sscanf(str,'%f,%f');
  p = NaN;
  t = NaN;
  if ( numel(vals) > 0 )
    p = vals(1);
  end;
  if ( numel(vals) > 1 )
    t = vals(2);
  end;
return;
