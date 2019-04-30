function review_avhrr_weekly(stn,dt,fld,deldt)
%function review_avhrr_weekly(stn,dt,fld,deldt)
%
% Plot a sequence of CONTOURF fields around station STN from STN.(FLD)
% (DEFAULT: 'raw_avhrr_weekly_sst_field'), for each of the DELDT (DEFAULT: 3)
% distinct dates prior to, closest to, and after datenum DT. (Note the normal
% gap between elements of STN.raw_avhrr_weekly_sst_field is *one week*.)
%
% Last Saved Time-stamp: <Thu 2012-03-29 15:12:09  Lew.Gramer>

  if ( ~exist('deldt','var') || isempty(deldt) )
      deldt=3;
  end;

  fld='raw_avhrr_weekly_sst_field';

  [dterr,dtix] = min(abs(stn.(fld).date-dt));
  if ( dterr > 14 )
      warning('Closest date displayed is >14 days from %s',datestr(dt));
  end;

  ixen = dtix-deldt:dtix+deldt;
  % Make color range consistent across all plots
  cmn = nanmin(nanmin(nanmin(stn.(fld).field(ixen,:,:))));
  cmx = nanmax(nanmax(nanmax(stn.(fld).field(ixen,:,:))));

  fh=fmg;
  for ix=ixen(:)'
      if ( ~ishandle(fh) )
        disp('Figure handle no longer valid...');
        break;
      end;
      figure(fh);
      hold off;
      contourf(stn.(fld).lon,stn.(fld).lat,squeeze(stn.(fld).field(ix,:,:)));
      hold on;
      caxis([cmn,cmx]);
      colorbar;
      plot(stn.lon,stn.lat,'kp');
      titlename([upper(stn.station_name),'.',strrep(fld,'_','\_'),' ',...
                 datestr(stn.(fld).date(ix)-3.5),' (',num2str(get_jday(stn.(fld).date(ix)-3.5)),')']);
      pause;
  end;
  if ( ishandle(fh) )
    close(fh);
  end;

return;
