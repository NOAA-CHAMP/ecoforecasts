function stns = read_NF09(stns)
%function stns = read_NF09(stns)
%
% Add the mooring data structs from Ecoforecasts/Data/'NF09_moorings.mat' as
% new fields in the struct STNS. If any field returned from the .MAT file is
% already in STNS, signal an error.) Assumes PLOT_NF09_MOORINGS_UPWELLING.m
% has previously been run, so that .MAT file exists.
%
% Last Saved Time-stamp: <Fri 2018-03-02 09:20:33 Eastern Standard Time gramer>

  if ( ~exist('stns','var') )
    stns = [];
  end;

  %nf09 = read_NF09_moorings;
  matfname = fullfile(get_ecoforecasts_path('data'),'NF09_moorings.mat');
  if ( exist(matfname,'file') )
    disp(['Loading ',matfname]);
    nf09 = load(matfname);
  else
    error('Missing %s!',matfname);
  end;

  flds = fieldnames(nf09);
  for fldix=1:numel(flds)
    fld = flds{fldix};
    if ( isfield(stns,fld) )
      error('STNS.%s already exists!');
    end;
    stns.(fld) = nf09.(fld);
  end;

return;
