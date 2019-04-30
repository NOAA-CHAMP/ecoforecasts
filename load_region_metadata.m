function stns = load_region_metadata(region)
%function stns = load_region_metadata(region)
%
% Load station codes, coordinates, and other metadata from the CSV file
% [REGION '_site_metadata.csv']. Returns a vector of structs STNS.
%
% Last Saved Time-stamp: <Fri 2011-04-01 10:34:17  Lew.Gramer>

  cfgpath = get_ecoforecasts_path;

  fname = fullfile(cfgpath,[region '_site_metadata.csv']);
  fid = fopen(fname, 'r');
  if ( fid < 0 )
    error('Unable to open metadata file "%s"!', fname);
  end;
  metad = textscan(fid,'%s %s %s %s %f %f %f','Delimiter',',');
  fclose(fid);

  commentix = find(strmatch('#',metad{1}));
  for ix = 1:size(metad,1)
    metad{ix}(commentix) = [];
  end;  

  for ix = 1:length(metad{1})
    stns(ix).station_name  = metad{1}{ix};
    stns(ix).name          = metad{2}{ix};
    stns(ix).fname         = metad{4}{ix};  % In situ data filename
    stns(ix).lon           = metad{5}(ix);
    stns(ix).lat           = metad{6}(ix);
    stns(ix).depth         = metad{7}(ix);
    stns(ix).misst_region  = region;
  end;

  % Basic QC
  for ix = length(stns):-1:1
    if ( any(isnan([stns(ix).lon stns(ix).lat stns(ix).depth])) )
      warning('Removing station record %d named "%s"...', ix, stns(ix).name);
      stns(ix) = [];
    end;
  end;

return;
