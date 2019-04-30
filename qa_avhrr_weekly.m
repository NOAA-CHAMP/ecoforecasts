1;
%%%%%%%%%% SCRIPT to do human QA on AVHRR Weekly SST at station STNM (DEFAULT: 'fwyf1')
%%%%%%%%%% 

datapath = get_ecoforecasts_path('data');

if ( ~exist('stnm','var') || isempty(stnm) )
  stnm = 'fwyf1';
end;

rad=0.5;

if ( ~exist('lms','var') || ~exist('stn','var') ...
     || ~isfield(stn,'lon') || ~isfield(stn,'lat') ...
     || ~isfield(stn,'station_name') || ~strcmpi(stn.station_name,stnm) )
  stn = []; clear stn
  stn = get_station_from_station_name(stnm);
  stn = load_all_ndbc_data(stn);
  stn = verify_variable(stn,'ndbc_sea_t_7_d_lp');

  %Pixel indices, e.g., for FWYF1 ~ [1150,1250,550,650], with set(gca,'ydir','rev'); 
  lms = [stn.lon-rad,stn.lon+rad,stn.lat-rad,stn.lat+rad];
end;

if ( ~isfield(stn,'avhrr_weekly_sst_field') || ~exist('fld','var') )
  fld=[]; clear fld
  stn = get_avhrr_weekly_field(stn,true);
  fld = stn.avhrr_weekly_sst_field;
end;

jds = 1:7:364;
ixen = find(ismember(get_jday(fld.date),jds));

ixen(get_year(fld.date(ixen))<1997) = [];

stn.avhrr_weekly_sst_max.date = fld.date(ixen);
stn.avhrr_weekly_sst_max.data = interp_field(fld.lat,fld.lon,fld.field(ixen,:,:),stn.lat,stn.lon,{@nanmax,3},nan);
stn.avhrr_weekly_sst_min.date = fld.date(ixen);
stn.avhrr_weekly_sst_min.data = interp_field(fld.lat,fld.lon,fld.field(ixen,:,:),stn.lat,stn.lon,{@nanmin,3},nan);

fh=fmg;
plot_ts(stn.ndbc_sea_t,'k-',stn.ndbc_sea_t_7_d_lp,':','Color',[.5,.5,.5],...
        stn.avhrr_weekly_sst_max,'rv',stn.avhrr_weekly_sst_min,'b^');

if ( ~exist('ok','var') )
  firstix = find(fld.date(ixen)>=datenum(1997,1,1),1);
  if ( isempty(firstix) )
    error('No dates in %s.avhrr_weekly_sst_field >= 1997-Jan-01??',upper(stnm));
  end;
  ok = repmat(0,[firstix-1,1]);
  % Assume first two weeks of 1997 are good
  ok(end+1) = 1;
  ok(end+1) = 1;
end;

for wk=1+2:numel(ixen)-2
  if ( ~ishandle(fh) )
    disp('Figure closed?');
    break;
  end;
  figure(fh);
  % xlim(fld.date(ixen([wk-2,wk+2])));
  xlim([fld.date(ixen(wk-2))-1,fld.date(ixen(wk+2))+1]);
  ylim('default');
  yl = mean(ylim);
  ylim([floor(yl)-2,ceil(yl)+2]);
  datetick3('x',2,'keeplimits');
  ok(end+1) = input(['Week #',num2str(wk),' (',datestr(fld.date(ixen(wk))),') OK? ']);
end;


yorn = input('Save QA flags? ','s');
if ( strcmpi(yorn(1),'y') )
  % Assume LAST two weeks of record are also good
  ok(end+1) = 1;
  ok(end+1) = 1;
  save(datapath,[lower(stnm),'_qa_avhrr_weekly.mat']);
end;




clear datapath jds ixen fh firstix wk yl yorn
return;




if ( ~exist('ok','var') )
  firstix = find(fld.date(ixen)>=datenum(1997,1,1),1);
  if ( isempty(firstix) )
    error('No dates in %s.avhrr_weekly_sst_field >= 1997-Jan-01??',upper(stnm));
  end;
  ok = repmat(0,[firstix-1,1]);
end;

if ( ~exist('LONS','var') )
  [ig,ig,ig,ig,ig,ig,LONS,LATS] = get_avhrr_coords('florida');
end;

for ix=numel(ok)+1:numel(ixen);
  ixix=ixen(ix);
  [jd,yr]=get_jday(fld.date(ixix));
  if ( jd == 358 )
    endjd = get_jday(datenum(yr,12,31));
  else
    endjd = jd+6;
  end;
  fname = sprintf('data/avhrr/all.%04d%03d.%04d%03d.sst.day_night.florida.anomaly.png',...
                  yr,jd,yr,endjd)
  fh=fmg;
  sst = imread(fname);
  sst(sst >= 254) = 0;
  imagesc(LONS,LATS,sst);
  plot(stn.lon,stn.lat,'kp');
  axis(lms);
  titlename(datestr(fld.date(ixix)));
  ok(end+1) = input(['OK ',num2str(ix),'? ']);
  close(fh);
end;
