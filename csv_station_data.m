function csv_station_data(stn,flds,fname,begdt,enddt,nonnan,forDistribution,hasBeenQCd)
%function csv_station_data(stn,flds,fname,begdt,enddt,nonnan,forDistribution)
%
% Save station structure data in "CSV" (comma-separated text) format. If FLDS
% is a nonempty string or cellstr, only save THOSE structure fields into the
% CSV file. If FNAME is a nonempty string, save CSV to that alternate path
% and filename. DEFAULT fname 'MATLABHOME/ecoforecasts/data/<stn>.csv'. If
% BEGDT and ENDDT are *both* valid datenums, then only data with timestamps
% between those two times (inclusive) are output to CSV file. If optional arg
% NONNAN==True, only output data for times when ALL output fields are non-NaN.
% If optional arg FORDISTRIBUTION is true, include title line, sensor units,
% and standard ICON Disclaimer as additional lines in the CSV file.
%
% Last Saved Time-stamp: <Tue 2018-03-06 11:13:07 EST lew.gramer>

  datapath = get_ecoforecasts_path('data');

  if ( ~exist('stn','var') || ~isfield(stn,'station_name') )
    error('First arg must be a station STRUCT with field "station_name"');
  end;
  if ( ~exist('flds','var') || isempty(flds) )
    flds = fieldnames(stn);
  end;
  if ( ischar(flds) )
    flds = {flds};
  end;

  if ( ~exist('fname','var') || isempty(fname) )
    fname = fullfile(datapath, [stn.station_name '.csv']);
  end;
  if ( ~ischar(fname) )
    error('Third arg FNAME must be name for CSV output file!');
  end;
  if ( exist(fname,'file') )
    warning('CSV file "%s" already exists: overwriting...', fname);
  end;

  if ( ~exist('begdt','var') || isempty(begdt) )
    begdt = -Inf;
  end;
  if ( ~exist('enddt','var') || isempty(enddt) )
    enddt = +Inf;
  end;

  if ( ~exist('nonnan','var') || isempty(nonnan) )
    nonnan = false;
  end;

  if ( ~exist('forDistribution','var') || isempty(forDistribution) )
    forDistribution = false;
  end;
  if ( ~exist('hasBeenQCd','var') || isempty(hasBeenQCd) )
    hasBeenQCd = false;
  end;


  dts = [];
  for ifld = 1:length(flds)
    fld = flds{ifld};
    if ( isfield(stn.(fld),'date') )
      dts = union(dts, stn.(fld).date);
    end;
  end;
  dts = dts(begdt <= dts & dts <= enddt);


  nRecords = 0;

  fid = fopen(fname, 'w');
  if ( fid < 0 )
    error('UNABLE TO OPEN CSV OUTPUT FILE "%s"!', fname);
  end;

  if ( forDistribution )
    if ( hasBeenQCd )
      fprintf(fid, '"NOAA CREWS/ICON Station %s"\n', stn.station_name);
      fprintf(fid, '"===================="\n');
      fprintf(fid, '"DISCLAIMER: These meteorological and oceanographic data have been screened for quality-control for accuracy, but are PRELIMINARY. NOAA can not be held liable for use of these data in any manner other than for the perusal of preliminary oceanographic data in scientific research on coral reefs. For ICON in situ data used in research or other publications, please credit the NOAA AOML Integrated Coral Observing Program (ICON/CREWS), Dr. J. C. Hendee, Principal Investigator. For all other data, please cite the appropriate partner institution as listed in our Sources web page for this station."\n');
      fprintf(fid, '"===================="\n');
      fprintf(fid, '\n');
      fprintf(fid, '\n');
    else
      fprintf(fid, '"NOAA CREWS/ICON Station %s"\n', stn.station_name);
      fprintf(fid, '"===================="\n');
      fprintf(fid, '"DISCLAIMER: These meteorological and oceanographic data are PRELIMINARY and have not been screened or quality-controlled for accuracy. NOAA can not be held liable for use of these data in any manner other than for the perusal of preliminary oceanographic data in scientific research on coral reefs. For ICON in situ data used in research or other publications, please credit the NOAA AOML Integrated Coral Observing Program (ICON/CREWS), Dr. J. C. Hendee, Principal Investigator. For all other data, please cite the appropriate partner institution as listed in our Sources web page for this station."\n');
      fprintf(fid, '"===================="\n');
      fprintf(fid, '\n');
      fprintf(fid, '\n');
    end;
  end;

  fprintf(fid, 'Station,Date,JYear,JDay,Hour,Minute,Second');
  for ifld = 1:length(flds)
    fld = flds{ifld};
    fprintf(fid, ',%s', fld);
  end;
  fprintf(fid, '\n');    

  for idt = 1:length(dts)
    dt = dts(idt);

    [yr,mo,dy,hr,mn,sc] = datevec(dt);
    jd = datenum(yr,mo,dy) - datenum(yr,1,1) + 1;

    dat = [];
    for ifld = 1:length(flds)
      fld = flds{ifld};
      if ( isfield(stn.(fld),'date') && isfield(stn.(fld),'data') )
        ix = find(stn.(fld).date == dt);
        if ( ~isempty(ix) )
          dat(ifld,1) = stn.(fld).data(ix);
        else
          dat(ifld,1) = -999.0;
        end;
      end;
    end;

    if ( ~nonnan || all(~isnan(dat(:))) )
      fprintf(fid, '%s,%s,%04d,%03d,%02d,%02d,%02d', stn.station_name, ...
              datestr(dt,0), yr, jd, hr, mn, round(sc));
      for ifld = 1:length(flds)
        fprintf(fid, ',%g', dat(ifld));
      end;
      fprintf(fid, '\n');
      nRecords = nRecords + 1;
    end;

  end;


  if ( forDistribution )
    fprintf(fid, '\n');
    fprintf(fid, '\n');
    fprintf(fid, '"===================="\n');
    fprintf(fid, 'UNITS FOR ICON DATA:\n');
    fprintf(fid, ',Wind speeds,Knots\n');
    fprintf(fid, ',Ocean current speeds,m s-1\n');
    fprintf(fid, ',Air pressure,hectoPascals\n');
    fprintf(fid, ',PAR (light 400-700nm),micro-mole quanta m-2 s-1\n');
    fprintf(fid, ',UV (light 305 or 330 or 380 nm),mW m-2 nm-1\n');
    fprintf(fid, ',Kd (attenuation coefficient),m-1\n');
    fprintf(fid, ',Temperatures (sea or air or dewpoint),deg-C\n');
    fprintf(fid, ',Seawater conductivity,milliSiemens cm-1\n');
    fprintf(fid, ',Sea pressure,decibars\n');
    fprintf(fid, ',Salinity,Practical salinity units (psu)\n');
    fprintf(fid, ',Water depth,meters\n');
    fprintf(fid, ',Sound velocity,m s-1\n');
    fprintf(fid, ',Chlorophyll a,micro-g m-2\n');
    fprintf(fid, ',Rain amount,mm hour-1\n');
    fprintf(fid, ',pCO2,micro-atmospheres\n');
    fprintf(fid, '"===================="\n');
    fprintf(fid, '\n');
  end;

  fclose(fid);

  disp(['Output ',num2str(nRecords),' records']);

return;
