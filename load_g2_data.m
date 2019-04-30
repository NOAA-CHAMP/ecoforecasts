function [dts, data] = load_g2_data(fname, begdt, enddt)
%function [dts, data] = load_g2_data(fname, begdt, enddt)
% Load data from a file in the format commonly used by the G2 TWNG Fancy
% Graphs 'export' feature. Assumes export was done comma-delimited. If the
% second and third args 'begdt' and 'enddt' specify a valid date range, only
% data that lies within that range is loaded from the file and returned.
% 
% Last Saved Time-stamp: <Wed 2008-08-27 09:42:05 Eastern Daylight Time gramer>

  dts = [];
  data = [];

  if ( ~exist('begdt') || isempty(begdt) )
    begdt = -inf;
  end;
  if ( ~exist('enddt') || isempty(enddt) )
    enddt = inf;
  end;

  fid = fopen(fname, 'r');

  if ( fid < 0 )
    error('Unable to open file %s', fname);
  end

  flds = { [], [], [], [], [], [] };

  % Parse each line of CSV(-ish) G2 export file

  if ( strfind(version('-release'), '200') )
    % Junk,Junk,Date Time,Value
    flds = textscan(fid, '%*[^,],%*f,%f/%f/%f %f:%f,%f');

  else
    % Version 6.x or earlier - have to do it the slooooow way
    while 1
      ln = fgetl(fid);
      if ( ~ischar(ln) ); break; end;
      % Junk,Junk,Date Time,Value
      idx = strfind(ln, ',');
      dtm = ln(idx(2)+1:idx(3)-1);
      flds{6} = [ flds{6} str2num(ln(idx(3)+1:end)) ];
      idx = strfind(dtm, '/');
      sidx = strfind(dtm, ' ');
      cidx = strfind(dtm, ':');
      flds{1} = [ flds{1} str2num(dtm(1:idx(1)-1)) ];
      flds{2} = [ flds{2} str2num(dtm(idx(1)+1:idx(2)-1)) ];
      flds{3} = [ flds{3} str2num(dtm(idx(2)+1:sidx-1)) ];
      flds{4} = [ flds{4} str2num(dtm(sidx+1:cidx-1)) ];
      flds{5} = [ flds{5} str2num(dtm(cidx+1:end)) ];
    end;
  end

  fclose(fid);


  rawdts = datenum(flds{3}, flds{1}, flds{2}, flds{4}, flds{5}, 0);
  rawdata = flds{6};

  [dts, uniqidx] = unique(rawdts);
  data = rawdata(uniqidx);

  data = data(begdt <= dts & dts <= enddt);
  dts = dts(begdt <= dts & dts <= enddt);

return;
