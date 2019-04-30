function station = extract_avhrr_weekly_field(station_or_stnm,begdt,enddt,dataset)
%function station = extract_avhrr_weekly_field(station_or_stnm,begdt,enddt,dataset)
%
% NOTE: Normally only called by GET_STATION_AVHRR_WEEKLY_FIELD, not directly.
%
% Extract a (default) 17x17 1km pixel field around the position of station
% contained in struct or named in string STATION_OR_STNM. Updates and returns
% station struct with new field STATION.AVHRR_WEEKLY_SST containing subfields
% scalar .LON and .LAT, Nx1 vector .DATE, 17x17 matrices .N and .PCTUSED, and
% Nx17x17 3-d matrix of SSTs STATION.AVHRR_WEEKLY_SST.FIELD. A once-a-day SST
% field is produced from mission inception (1993) to the present, using AVHRR
% Weekly composites with a INTERP1 shape-preserving spline ('pchip').
%
% Last Saved Time-stamp: <Tue 2012-07-03 14:35:30  lew.gramer>

  %DEBUG:
  tic,
  station = get_station_from_station_name(station_or_stnm);
  station_or_stnm = []; clear station_or_stnm;

  % all.1993239.1993245.sst.day_night.florida.mean.png has heavy cloudcover
  zerodate = datenum(1993,1,1) + 246 - 1;

  if ( ~exist('begdt','var') || isempty(begdt) )
    begdt = zerodate;
  end;
  if ( ~exist('enddt','var') || isempty(enddt) )
    enddt = floor(now) - 14;
  end;
  if ( ~exist('dataset','var') || isempty(dataset) )
    dataset = 'mean';
  end;

  % This is the way USF names/numbers their weekly composite PNG files
  wkjds = 1:7:364;

  [begyr,begmo,begdy] = datevec(begdt);
  [endyr,endmo,enddy] = datevec(enddt);

  dts = [];
  for yr = begyr:endyr
    dts(end+1:end+length(wkjds)) = datenum(yr,1,1) + wkjds - 1;
  end;
  dts(dts < (begdt-6)) = [];
  dts(dts > (enddt+7)) = [];


  [ig,LON,LAT] = query_avhrr_weekly_subset(station,'bounds');

  rowrad = 8;
  colrad = 8;
  midrow = round(size(LON,1)/2);
  midcol = round(size(LON,2)/2);
  rows = midrow-rowrad:midrow+rowrad;
  cols = midcol-colrad:midcol+colrad;

  station.avhrr_weekly_sst.lon = LON(1,cols)';
  station.avhrr_weekly_sst.lat = LAT(rows,1);

  sst = repmat( nan, [length(dts) length(rows) length(cols)] );

  for dtix = 1:length(dts)
    dt = dts(dtix);
    [yr,mo,dy] = datevec(dt);
    jd = dt - datenum(yr,1,1) + 1;
    wk = ceil(jd/7);
    wholesst = query_avhrr_weekly_subset(station,yr,wk,dataset);
    sst(dtix,:,:) = wholesst(rows,cols);
  end;

  %DEBUG:
  toc,

  %DEBUG:  tic,
  alldts = dts(1):dts(end);
  station.avhrr_weekly_sst.date = alldts';
  station.avhrr_weekly_sst.N = repmat(nan, [size(sst,2) size(sst,3)]);
  station.avhrr_weekly_sst.pctused = repmat(nan, [size(sst,2) size(sst,3)]);
  station.avhrr_weekly_sst.field = repmat(nan, [numel(alldts) size(sst,2) size(sst,3)]);

  % Quick and dirty way of finding all non-land pixels (ASSUMING sst has
  % enough samples to be a reliable indicator of land vs. sea!)
  meansst = squeeze(nanmean(sst,1));
  seaix = find(isfinite(meansst));

  nwks = size(sst,1);
  for ix = seaix(:)'
    usedix = find(isfinite(sst(:,ix)));
    station.avhrr_weekly_sst.N(ix) = numel(usedix);
    station.avhrr_weekly_sst.pctused(ix) = (numel(usedix)/nwks);
    station.avhrr_weekly_sst.field(:,ix) = ...
        interp1(dts(usedix),sst(usedix,ix),alldts,'pchip',nan);
  end;
  %DEBUG:  toc,

  sst = []; clear sst;

return;
