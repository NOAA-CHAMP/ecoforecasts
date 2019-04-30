function stn = load_gfs_station(stnm,mdl,grd)
%function stn = load_gfs_station(stnm[,mdl[,grd]])

    stnm = lower(stnm);
    if ( ~exist('mdl','var') || isempty(mdl) )
        mdl = 'gfs';
    end;
    if ( ~exist('grd','var') || isempty(grd) )
        grd = '0p25';
    end;
    stn = [];
    
    mdlgrd = [mdl,'_',grd];
    
    thisyear = get_year(now);
    for yr=2015:thisyear
        matfname = [mdlgrd,'_',stnm,'_',num2str(yr),'.mat'];
        load(matfname);
        fldnms = fieldnames(stns.(stnm));
        for fldix=1:numel(fldnms)
            fldnm = fldnms{fldix};
            sfldnm = [mdlgrd,'_',fldnm];
            if ( is_ts(stns.(stnm).(fldnm)) )
                stns.(stnm).(sfldnm).date = stns.(stnm).(fldnm).date(:);
                stns.(stnm).(sfldnm).data = stns.(stnm).(fldnm).data(:);
                stns.(stnm) = rmfield(stns.(stnm),fldnm);
            end;
        end;
        stn = merge_station_data(stn,stns.(stnm));
    end;

return;
