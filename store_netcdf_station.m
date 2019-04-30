function store_netcdf_station(parid,stn,stnm,attrs,longstnm,vars)
%function store_netcdf_station(parid,stn,stnm,attrs,longstnm,vars)
%
% Store a netCDF Group ID for the station STN with name STNM in the parent
% (or global) netCDF Group PARID. Each variable in VARS is stored in the new
% group as its own sub-group.
%
% Last Saved Time-stamp: <Sun 2018-03-04 17:58:46 Eastern Standard Time gramer>

  if ( isempty(stnm) && isfield(stn,'station_name') )
    stnm = stn.station_name;
  end;

  NC_GLOBAL = netcdf.getConstant('NC_GLOBAL');
  fillValue = -9999;

  % netcdf.reDef(ncid);

  gid = netcdf.defGrp(parid,lower(stnm));
  netcdf.putAtt(gid,NC_GLOBAL,'name',upper(stnm));
  netcdf.putAtt(gid,NC_GLOBAL,'long_name',longstnm);

  if ( isfield(stn,'lon') )
    netcdf.putAtt(gid,NC_GLOBAL,'longitude',stn.lon);
  end;
  if ( isfield(stn,'lat') )
    netcdf.putAtt(gid,NC_GLOBAL,'latitude',stn.lat);
  end;
  if ( isfield(stn,'depth') )
    netcdf.putAtt(gid,NC_GLOBAL,'depth',stn.depth);
  end;
  if ( isfield(stn,'isobath_orientation') )
    netcdf.putAtt(gid,NC_GLOBAL,'isobath_angle',stn.isobath_orientation);
  end;
  if ( isfield(stn,'slope') )
    netcdf.putAtt(gid,NC_GLOBAL,'seafloor_slope',stn.slope);
  end;

  for attix = 1:numel(attrs)
    attr = attrs{attix};
    if ( isfield(stn,attr) )
      netcdf.putAtt(gid,NC_GLOBAL,attr,stn.(attr));
    end;
  end;


  for varix = 1:numel(vars)
    dps = [];
    var = vars{varix};
    desc = '';

    ts = [];
    if ( iscell(var) )
      if ( numel(var) < 2 || numel(var) > 3 || ~ischar(var{2}) )
        error(['If arg VARS is a cell, it may have optional numeric depth vector *or* time series STRUCT, then CHAR field name, and then an optional description string']);
      end;
      if ( isnumeric(var{1}) )
        dps = var{1};
        var(1) = [];
      elseif ( numel(var) == 3 && ischar(var{2}) )
        if ( ~is_valid_ts(var{1}) )
          error('Invalid time series STRUCT passed in for STN.(%s)',var{2});
        end;
        ts = var{1};
        var(1) = [];
      end;
      if ( numel(var) > 1 )
        if ( ~ischar(var{2}) )
          error('Optional description cell elt. must be a CHAR');
        else
          desc = var{2};
        end;
      end;
      var = var{1};
    end;
    if ( ~ischar(var) )
      error('VAR arg #%d must specify a CHAR fieldname in %s',varix,stnm);
    end;
    if ( isempty(ts) )
      ts = stn.(var);
    end;

    var_gid = netcdf.defGrp(gid,lower(var));
    netcdf.putAtt(var_gid,NC_GLOBAL,'Description',desc);
    var_tim_dim_id = netcdf.defDim(var_gid,'time',[numel(ts.date)]);
    var_dts_id = netcdf.defVar(var_gid,'date','double',[var_tim_dim_id]);
    if ( isfield(ts,'data') )
      var_dat_id = netcdf.defVar(var_gid,'data','double',[var_tim_dim_id]);
    end;
    var_dep_id = [];
    if ( isfield(ts,'prof') )
      if ( size(ts.prof,2) ~= numel(dps) )
        warning('Ignoring %s.%s.prof field: no valid depth vector specified',stnm,var);
      else
        var_dep_dim_id = netcdf.defDim(var_gid,'depth',[numel(dps)]);
        var_dep_id = netcdf.defVar(var_gid,'depth','double',[var_dep_dim_id]);
        var_prf_id = netcdf.defVar(var_gid,'prof','double',[var_dep_dim_id,var_tim_dim_id]);
      end;
    end;

    % netcdf.unDef(ncid);
    netcdf.putVar(var_gid,var_dts_id,ts.date);
    if ( isfield(ts,'data') )
      netcdf.putVar(var_gid,var_dat_id,ts.data);
    end;
    if ( ~isempty(var_dep_id) )
      netcdf.putVar(var_gid,var_dep_id,dps);
      netcdf.putVar(var_gid,var_prf_id,ts.prof);
    end;
    % netcdf.reDef(ncid);
  end;

return;
