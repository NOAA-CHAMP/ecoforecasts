function res = get_sfomc_trim_data
%function res = get_sfomc_trim_data
%
% Get data from the South Florida Oil Monitoring Center field project of Nova
% Southeastern University's Oceanography Center (Soloviev et al.) NOTE: This
% version REMOVES several memory-intensive fields that are rarely needed.
%
% Last Saved Time-stamp: <Tue 2016-08-30 16:46:51 Eastern Daylight Time lew.gramer>

  datapath = get_ecoforecasts_path('data');

  trim_matfname = fullfile(datapath,'SFOMC_trim.mat');

  if ( exist(trim_matfname,'file') )
    disp(['Loading ',trim_matfname]);
    load(trim_matfname);

  else
    res = get_sfomc_data;


    disp('Trimming SFOMC');

    % % Error determinants
    % res.nw_w_btm = rmfield(res.nw_w_btm,{'adcp_backscatter','adcp_cacount','adcp_err'});

    % Miscellaneous fields
    res.nw_w_btm = rmfield(res.nw_w_btm,grepstruct(res.nw_w_btm,'_h_'));
    res.ne_buoy = rmfield(res.ne_buoy,grepstruct(res.ne_buoy,'delta_seatemp'));

    % These can be back-calculated: Raw direction is hard to work with anyway
    res.nw_w_btm = rmfield(res.nw_w_btm,grepstruct(res.nw_w_btm,'_([uv]|dir)$'));
    res.nw_w_btm = rmfield(res.nw_w_btm,grepstruct(res.nw_w_btm,'_([uv]|dir)_'));
    res.ne_buoy = rmfield(res.ne_buoy,grepstruct(res.ne_buoy,'_([uv]|dir)$'));
    res.ne_buoy = rmfield(res.ne_buoy,grepstruct(res.ne_buoy,'_([uv]|dir)_'));

    % % Baroclinic profiles and time series
    % res.nw_w_btm = rmfield(res.nw_w_btm,grepstruct(res.nw_w_btm,'baroclinic'));
    % res.ne_buoy = rmfield(res.ne_buoy,grepstruct(res.ne_buoy,'baroclinic'));

    % % Dump baroclinic profiles, but keep depth-averaged baroclinic time series
    % res.nw_w_btm = rmfield(res.nw_w_btm,grepstruct(res.nw_w_btm,'baroclinic_[^_]$'));
    % res.ne_buoy = rmfield(res.ne_buoy,grepstruct(res.ne_buoy,'baroclinic_[^_]$'));

    % % Sub-sample baroclinic fields
    % flds = grepstruct(res.nw_w_btm,'_baroclinic');

    % Sub-sample all ADCP time series fields to hourly values
    flds = grepstruct(res.nw_w_btm,'adcp_');
    for fldix = 1:numel(flds)
      fld = flds{fldix};
      ts1 = res.nw_w_btm.(fld);
      if ( isfield(ts1,'date') )
        newts.date = unique(round(ts1.date*24)/24);
        if ( isfield(ts1,'data') )
          newts.data = interp1(ts1.date,ts1.data,newts.date,'nearest');
          if ( numel(newts.data) ~= numel(newts.date) )
            error('Failed to sub-sample RES.NW_W_BTM.(%s).data',fld);
          end;
        end;
        if ( isfield(ts1,'prof') )
          newts.prof = interp1(ts1.date,ts1.prof,newts.date,'nearest');
          if ( size(newts.prof,1) ~= numel(newts.date) )
            error('Failed to sub-sample RES.NW_W_BTM.(%s).prof',fld);
          end;
        end;
        res.nw_w_btm = rmfield(res.nw_w_btm,fld);
        res.nw_w_btm.(fld) = newts;
        newts=[]; clear newts
      end;
      ts1=[]; clear ts1
    end;
    clear fldix flds fld

    disp(['Saving ',trim_matfname]);
    save(trim_matfname,'res','-v7.3');

  end; %if ( exist(trim_matfname,'file') ) else

return;
