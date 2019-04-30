1;

error('RUN-ONCE HACK SCRIPT! Do not run this...');

  clear fwyf1; clear mlrf1;

  fwyf1 = load_all_ndbc_data([],'fwyf1');
  mlrf1 = load_all_ndbc_data([],'mlrf1');

  fwyf1 = verify_variable(fwyf1,'ndbc_wind1_u');
  fwyf1 = verify_variable(fwyf1,'ndbc_wind1_v');
  mlrf1 = verify_variable(mlrf1,'ndbc_wind1_u');
  mlrf1 = verify_variable(mlrf1,'ndbc_wind1_v');

  flds = grepstruct(fwyf1,'ndbc_');

  for cstnm=grepstruct(stns,'bnp')';
    stnm=cstnm{:};
    disp(stnm);
    station.station_name=upper(stnm);
    % station=stns.(stnm);
    % station=rmfield(station,'seatemp');

    for cfld=flds'
      fld=cfld{:};
      if ( isfield(mlrf1,fld) )
        disp(fld);
        %% HACK??? FWYF1 and MLRF1 instruments are at different HEIGHTS!
        % station.(fld) = ts_op(fwyf1.(fld),mlrf1.(fld),@(x,y)(((2.*x)+(y))./3));
        station.(fld) = fwyf1.(fld);
      end;
    end;

    station.ndbc_wind1_dir.date = station.ndbc_wind1_u.date;
    station.ndbc_wind1_dir.data = uv_to_dir(station.ndbc_wind1_u.data,station.ndbc_wind1_v.data);

    fname=fullfile('data',[lower(station.station_name),'-ndbc.mat']),
    save(fname,'station');
    station=[]; clear station;
  end;

  clear fwyf1 mlrf1 cfld cstnm fld flds fname stnm
