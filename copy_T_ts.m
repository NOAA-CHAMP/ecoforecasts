function ts = copy_T_ts(stns,ts,lbl)
%function ts = copy_T_ts(stns,ts,lbl)

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

    stix = numel(ts) + 1;
    ts(stix).station_name = [lbl,upper(stnm)];
    ts(stix).lon = stns.(stnm).lon;
    ts(stix).lat = stns.(stnm).lat;
    ts(stix).depth = stns.(stnm).depth;
    if ( isfield(stns.(stnm),'ngdc_depth') )
      ts(stix).ngdc_depth = stns.(stnm).ngdc_depth;
    end;
    if ( isfield(stns.(stnm),'ngdc_beta') )
      ts(stix).ngdc_beta = stns.(stnm).ngdc_beta;
    elseif ( isfield(stns.(stnm),'slope') )
      ts(stix).ngdc_beta = stns.(stnm).slope;
    else
      keyboard;
    end;
    if ( isfield(stns.(stnm),'ngdc_beta_ang') )
      ts(stix).ngdc_beta_ang = stns.(stnm).ngdc_beta_ang;
    elseif ( isfield(stns.(stnm),'slope_orientation') )
      ts(stix).ngdc_beta_ang = stns.(stnm).slope_orientation;
    else
      keyboard;
    end;
    if ( isfield(stns.(stnm),'seatemp') )
      ts(stix).date = stns.(stnm).seatemp.date;
      ts(stix).data = stns.(stnm).seatemp.data;
      if ( isfield(stns.(stnm).seatemp,'prof') )
        ts(stix).prof = stns.(stnm).seatemp.prof;
      end;
    elseif ( isfield(stns.(stnm),'t') )
      ts(stix).date = stns.(stnm).t.date;
      ts(stix).data = stns.(stnm).t.data;
    elseif ( isfield(stns.(stnm),'adcp_seatemp') )
      ts(stix).date = stns.(stnm).adcp_seatemp.date;
      ts(stix).data = stns.(stnm).adcp_seatemp.data;
    else
      keyboard;
    end;
  end;

return;
