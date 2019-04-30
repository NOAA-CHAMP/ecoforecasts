function sites = copy_T_sites(stns,sites,lbl)
%function sites = copy_T_sites(stns,sites,lbl)

  if ( iscell(stns) )
    cstns = stns; stns = [];
    for stnix=1:numel(cstns);
      stnm = genvarname(cstns{stnix}.station_name);
      stns.(stnm) = cstns{stnix};
    end;
    cstns=[]; clear cstns
  end;
  if ( ~exist('lbl','var') )
    lbl = '';
  end;

  %for str={'face','sefcri','sfomc'};
  cstnms = fieldnames(stns);
  for cstnmix = 1:numel(cstnms);
    stnm = cstnms{cstnmix};
    if ( ~isfield(stns.(stnm),'lon') )
      continue;
    end;

    stix = numel(sites) + 1;
    sites(stix).station_name = [lbl,upper(stnm)];
    sites(stix).lon = stns.(stnm).lon;
    sites(stix).lat = stns.(stnm).lat;
    if ( isfield(stns.(stnm),'depth') )
      sites(stix).depth = stns.(stnm).depth;
    else
      warning('No "%s" DEPTH',stnm);
    end;
    if ( isfield(stns.(stnm),'ngdc_depth') )
      sites(stix).ngdc_depth = stns.(stnm).ngdc_depth;
    else
      warning('No "%s" NGDC_DEPTH',stnm);
    end;
    if ( isfield(stns.(stnm),'ngdc_beta') )
      sites(stix).ngdc_beta = stns.(stnm).ngdc_beta;
    elseif ( isfield(stns.(stnm),'slope') )
      sites(stix).ngdc_beta = stns.(stnm).slope;
    else
      warning('No "%s" NGDC_BETA',stnm);
    end;
    if ( isfield(stns.(stnm),'ngdc_beta_ang') )
      sites(stix).ngdc_beta_ang = stns.(stnm).ngdc_beta_ang;
    elseif ( isfield(stns.(stnm),'slope_orientation') )
      sites(stix).ngdc_beta_ang = stns.(stnm).slope_orientation;
    else
      warning('No "%s" NGDC_BETA_ANG',stnm);
    end;

    % Sea temperature
    if ( isfield(stns.(stnm),'seatemp') )
      sites(stix).seatemp_fldnm = 'seatemp';
      sites(stix).seatemp = stns.(stnm).seatemp;
    elseif ( isfield(stns.(stnm),'t') )
      sites(stix).seatemp_fldnm = 't';
      sites(stix).seatemp = stns.(stnm).t;
    elseif ( isfield(stns.(stnm),'adcp_seatemp') )
      sites(stix).seatemp_fldnm = 'adcp_seatemp';
      sites(stix).seatemp = stns.(stnm).adcp_seatemp;
    else
      warning('No "%s" SEATEMP',stnm);
    end;

    % Meteorology
    flds = fieldnames(stns.(stnm));
    %windix = regexp(flds,'wind');
    %if ( any(~cellfun(@isempty,windix)) )
    windix = find(~cellfun(@isempty,regexp(flds,'wind')));
    if ( ~isempty(windix) )
      %spdix = regexp(flds(windix),'(speed|spd)');
      spdix = find(~cellfun(@isempty,regexp(flds(windix),'(speed|spd)')));
      if ( ~isempty(spdix) )
        spdix = windix(spdix);
        fld = flds{spdix};
        sites(stix).windspd_fldnm = fld;
        sites(stix).windspd = stns.(stnm).(fld);
      else
        warning('No "%s" WINDSPD',stnm);
        %%%%keyboard;
      end;
      %dirix = regexp(flds(windix),'dir');
      dirix = find(~cellfun(@isempty,regexp(flds(windix),'dir')));
      if ( ~isempty(dirix) )
        dirix = windix(dirix);
        fld = flds{dirix};
        sites(stix).winddir_fldnm = fld;
        sites(stix).winddir = stns.(stnm).(fld);
      else
        warning('No "%s" WINDDIR',stnm);
        %%%%keyboard;
      end;
    end;

    % airix = regexp(flds,'air');
    % if ( any(~cellfun(@isempty,airix)) )
    airix = find(~cellfun(@isempty,regexp(flds,'air')));
    if ( ~isempty(airix) )
      fld = flds{airix};
      sites(stix).airtemp_fldnm = fld;
      sites(stix).airtemp = stns.(stnm).(fld);
    elseif ( ~isempty(windix) )
      % IGNORE
    else
      warning('No "%s" AIRTEMP',stnm);
      %%%%keyboard;
    end;

  end;

return;
