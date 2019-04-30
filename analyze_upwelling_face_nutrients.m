face_path = get_coral_path('FACE');

for ccN = {'NF08','NF09'};
  cN = ccN{:};
  disp(cN);

  if ( ~exist('face','var') || ~isfield(face,cN) )
    face.(cN).cruiseName = cN;
  end;

  if ( strcmpi(cN,'NF08') )

    face.(cN).fname = fullfile(face_path,'Master_Data_Sheet_21aug09.xlsx');
    xl = importdata(face.(cN).fname);

    face.(cN).dt = [];
    face.(cN).stnm = {};
    face.(cN).lat = [];
    face.(cN).lon = [];
    face.(cN).z = [];
    face.(cN).t = [];
    face.(cN).s = [];
    face.(cN).tss = [];
    face.(cN).chl = [];
    face.(cN).pha = [];
    face.(cN).n = [];
    face.(cN).p = [];
    face.(cN).si = [];
    for csht = {'EddyXmptE1','EddyXmptE2'};
      sht = csht{:};
      disp(sht);
      endoff = 5;
      N = length(xl.data.(sht)(3:end-endoff,1));
      face.(cN).dt(end+1:end+N,1)   = datenum(xl.textdata.(sht)(3:end-endoff,1));
      face.(cN).stnm(end+1:end+N,1) = xl.textdata.(sht)(3:end-endoff,2);
      face.(cN).lat(end+1:end+N,1)  = xl.data.(sht)(3:end-endoff,1);
      face.(cN).lon(end+1:end+N,1)  = xl.data.(sht)(3:end-endoff,2);
      face.(cN).z(end+1:end+N,1)    = -xl.data.(sht)(3:end-endoff,3);
      face.(cN).t(end+1:end+N,1)    = xl.data.(sht)(3:end-endoff,4);
      face.(cN).s(end+1:end+N,1)    = xl.data.(sht)(3:end-endoff,5);
      face.(cN).tss(end+1:end+N,1)  = xl.data.(sht)(3:end-endoff,6);
      face.(cN).chl(end+1:end+N,1)  = xl.data.(sht)(3:end-endoff,7);
      face.(cN).pha(end+1:end+N,1)  = xl.data.(sht)(3:end-endoff,8);
      face.(cN).n(end+1:end+N,1)    = xl.data.(sht)(3:end-endoff,9);
      face.(cN).p(end+1:end+N,1)    = xl.data.(sht)(3:end-endoff,10);
      face.(cN).si(end+1:end+N,1)   = xl.data.(sht)(3:end-endoff,11);
    end;

    badix = find(isnan(face.(cN).z) | isnan(face.(cN).t));
    face.(cN).dt(badix) = [];
    face.(cN).stnm(badix) = [];
    face.(cN).lat(badix) = [];
    face.(cN).lon(badix) = [];
    face.(cN).z(badix) = [];
    face.(cN).t(badix) = [];
    face.(cN).s(badix) = [];
    face.(cN).tss(badix) = [];
    face.(cN).chl(badix) = [];
    face.(cN).pha(badix) = [];
    face.(cN).n(badix) = [];
    face.(cN).p(badix) = [];
    face.(cN).si(badix) = [];
    clear badix

  elseif ( strcmpi(cN,'NF09') )

    face.(cN).fname = fullfile(face_path,'Master_Data_Sheet_01Feb10.xls');
    xl = importdata(face.(cN).fname);

    sht = 'Sheet3';
    endoff = 0; dtendoff = endoff+21;
    face.(cN).dt = datenum(xl.textdata.(sht)(4:end-dtendoff,1))+xl.data.(sht)(1:end-endoff,1);
    face.(cN).stnm = xl.textdata.(sht)(4:end-dtendoff,3);
    face.(cN).lat = xl.data.(sht)(1:end-endoff,3)+(xl.data.(sht)(1:end-endoff,4)./60);
    face.(cN).lon = -( xl.data.(sht)(1:end-endoff,5)+(xl.data.(sht)(1:end-endoff,6)./60) );
    face.(cN).z = -xl.data.(sht)(1:end-endoff,7);
    face.(cN).t = xl.data.(sht)(1:end-endoff,8);
    face.(cN).s = xl.data.(sht)(1:end-endoff,9);
    face.(cN).chl = xl.data.(sht)(1:end-endoff,11);
    face.(cN).pha = xl.data.(sht)(1:end-endoff,12);
    face.(cN).n = xl.data.(sht)(1:end-endoff,13);
    face.(cN).no2 = xl.data.(sht)(1:end-endoff,14);
    face.(cN).no3 = face.(cN).n-face.(cN).no2;
    face.(cN).nh4 = xl.data.(sht)(1:end-endoff,15);
    face.(cN).p = xl.data.(sht)(1:end-endoff,16);
    face.(cN).si = xl.data.(sht)(1:end-endoff,17);

    badix = find(isnan(face.(cN).z) | isnan(face.(cN).t));
    face.(cN).dt(badix) = [];
    face.(cN).stnm(badix) = [];
    face.(cN).lat(badix) = [];
    face.(cN).lon(badix) = [];
    face.(cN).z(badix) = [];
    face.(cN).t(badix) = [];
    face.(cN).s(badix) = [];
    face.(cN).chl(badix) = [];
    face.(cN).pha(badix) = [];
    face.(cN).n(badix) = [];
    face.(cN).no2(badix) = [];
    face.(cN).no3(badix) = [];
    face.(cN).nh4(badix) = [];
    face.(cN).p(badix) = [];
    face.(cN).si(badix) = [];
    clear badix

  else
    error('Which cruise??');
  end;

  face.(cN).den = sw_dens(face.(cN).s,face.(cN).t,-face.(cN).z);
  face.(cN).dena = face.(cN).den - 1000;

  xl=[]; clear xl sht endoff dtendoff

  face.(cN).zmin   = -290;             face.(cN).zmax   = 0;
  face.(cN).tmin   = 7.9;              face.(cN).tmax   = 29.1;
  face.(cN).smin   = 34.9;             face.(cN).smax   = 36.8;
  face.(cN).nmin   = 0;                face.(cN).nmax   = 40;
  face.(cN).no2min = 0;                face.(cN).no2max = 1.6;
  face.(cN).nh4min = 0;                face.(cN).nh4max = 20;
  face.(cN).chlmin = 0;                face.(cN).chlmax = 0.65;
  face.(cN).phamin = 0;                face.(cN).phamax = 0.35;
  face.(cN).pmin   = 0;                face.(cN).pmax   = 2.0;
  face.(cN).simin  = 0;                face.(cN).simax  = 20.0;
  face.(cN).dmin   = 0;                face.(cN).dmax   = 27;

  %if ( ~exist('h','var') )
  if ( ~isfield(face.(cN),'h') )
    coastpath = get_ecoforecasts_path('coast');
    nc = mDataset(fullfile(coastpath,'fl_east_gom_crm_v1.nc'));
    if ( ~isempty(nc) )
      hlon = cast(nc{'x'}(:),'double');
      hlat = cast(nc{'y'}(:),'double');
      %xix = find(min(face.(cN).lon)-0.1<=hlon & hlon<=max(face.(cN).lon)+0.1);
      %yix = find(min(face.(cN).lat)-0.1<=hlat & hlat<=max(face.(cN).lat)+0.1);
      %hz = cast(nc{'z'}(yix(1):yix(end),xix(1):xix(end)),'double');
      zvar = nc{'z'};
      face.(cN).h = repmat(nan,size(face.(cN).lon));
      for dix = 1:numel(face.(cN).lon)
        [lonerr,lonix] = min(abs(hlon-face.(cN).lon(dix)));
        [laterr,latix] = min(abs(hlat-face.(cN).lat(dix)));
        if ( lonerr > (190/111e3) || laterr > (190/111e3) )
          clear zvar; close(nc); clear nc;
          error('Could not locate data point %d in NGDC 3" data: %g,%g',hix,face.(cN).lon,face.(cN).lat);
        end;
        face.(cN).h(dix) = cast(zvar(latix,lonix),'double');
      end;
      zvar = []; clear zvar; close(nc);
    end;
    clear nc;
    hlat = [];
    hlon = [];
    clear coastpath dix hlat hlon laterr latix lonerr lonix
  end; %if ( ~isfield(face.(cN),'h') )

  % This should include ALL "eddy" stations from NF08
  face.(cN).nfrt_ix = find(face.(cN).lat>25);


  % NOTE: Using "eddy" stations from NF08 already filters out boil stations, so
  % we have no need to filter based on station name for that cruise...

  %non_boil_stnix = union(strmatch('PE',face.(cN).stnm),strmatch('B',face.(cN).stnm));
  face.(cN).non_boil_stnix = union(strmatch('PE',face.(cN).stnm),strmatch('B',face.(cN).stnm));


  if ( ~isempty(face.(cN).non_boil_stnix) )
    face.(cN).below_25_ix = find(face.(cN).t(face.(cN).non_boil_stnix)<25);
    face.(cN).below_25_ix = face.(cN).non_boil_stnix(face.(cN).below_25_ix);
  else
    face.(cN).below_25_ix = find(face.(cN).t<25);
  end;

  if ( ~isempty(face.(cN).non_boil_stnix) )
    face.(cN).below_24_ix = find(face.(cN).t(face.(cN).non_boil_stnix)<=24);
    face.(cN).below_24_ix = face.(cN).non_boil_stnix(face.(cN).below_24_ix);
  else
    face.(cN).below_24_ix = find(face.(cN).t<=24);
  end;

  clear ans csht N;

end; %for ccN = {'NF08','NF09'};
