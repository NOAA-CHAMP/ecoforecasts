function stn = read_station_ef_csv(stn_or_stnm,efnm,efsub,fname)
%function stn = read_station_ef_csv(stn_or_stnm,efnm,efsub,fname)
%
% Read CSV (text) file exported from G2 for station named STN.station_name or
% STNM, and containing details of Ecoforecast alerts (G2 "imn-events") from
% an EF model named EFNM. Optionally, only load alerts for the EF *sub-model*
% named EFSUB. Return as a set of time series fields (date and S/RI) in STN.
% DEFAULT: Reads file GET_ECOFORECASTS_PATH('data'),[STNM,'_',EFNM,'.csv'].
% Alternate full or relative pathname for CSV file may be given in FNAME.
%
% Last Saved Time-stamp: <Tue 2013-07-23 15:29:27 Eastern Daylight Time gramer>

    stn = get_station_from_station_name(stn_or_stnm);
    if ( ~exist('fname','var') || isempty(fname) )
      fname = fullfile(get_ecoforecasts_path('data'),[upper(stn.station_name),'_',upper(efnm),'.csv']);
    end;
    [dt,yr,jd,st,typ,cls,subtyp,SRI] = textread(fname,'%[^,],%d,%d,%[^,],%[^,],%[^,],%[^,],%d,%*[^\n]\n','headerlines',1);
    dts = datenum(dt,'yyyy-mm-dd');

    if ( exist('efsub','var') )
        efnm = efsub;
        efix = find(strcmpi(efsub,subtyp));
        dts=dts(efix);
        SRI=SRI(efix);
    end;

    varnm = strrep(efnm,'-','_');

    udts=unique(dts);
    stn.(varnm).sri.date = udts(:);
    for ix=1:numel(udts);
        ixix=find(dts==udts(ix));
        if (isempty(ixix))
            stn.(varnm).sri.data(ix,1)=nan;
        else
            stn.(varnm).sri.data(ix,1)=sum(SRI(ixix));
        end;
    end;

return;
