function [annual, amj, jas, ond, jfm] = seasonal_data_ranges(stn, fieldname, validrange, validhours)
%FUNCTION [annual, amj, jas, ond, jfm] = seasonal_data_ranges(stn, fieldname, validrange, validhours)
%
% Categorize a time series into value ranges appropriate for each of four seasons...
%
% INPUTS:
%  stn - struct with station data
%  fieldname - name of a field in "stn" containing data to be categorized
%  validrange - [ minimum-acceptable-value maximum-acceptable-value ]
%      DEFAULT: All property values are acceptable
%  validhours - [ earliest-hour-of-day-to-consider latest-hour-to-consider ]
%      DEFAULT: All hours are acceptable
%
% OUTPUTS:
%  See 'help derive_prop_ranges' for the layout of each of the five return
%  structures - annual, and one each for Apr-May-Jun (amj), Jul-Aug-Sep (jas),
%  Oct-Nov-Dec (ond) and Jan-Feb-Mar (jfm), resp.
%
% DEPENDENCIES:
%  derive_prop_ranges.m - CLIPS-specific statistical treatment of a dataset.
%
% Last Saved Time-stamp: <Thu 2008-09-11 12:27:04 Eastern Daylight Time gramer>

    if ( ~isfield(stn, fieldname) )
        error('"%s" is not a field of "stn"!', fieldname);
    end
    if ( ~isfield(stn.(fieldname), 'date') ...
         || ~isfield(stn.(fieldname), 'data') )
        error('"stn.%s" needs both a date and a data field!', fieldname);
    end
    if ( ~exist('validrange', 'var') )
        validrange = [-Inf +Inf];
    end
    if ( ~exist('validhours', 'var') )
        validhours = [0 24];
    end

    annual = [];
    amj = [];
    jas = [];
    ond = [];
    jfm = [];


    dts = stn.(fieldname).date;
    dat = stn.(fieldname).data;

    [Y M D h m s] = datevec(dts);
    jday = floor(dts - datenum(Y, 1, 1) + 1);

    % Extract only those rows from the times of day we want to consider...
    goodhours = (validhours(1) <= h & h <= validhours(2));
    % And extract only those rows which contain a "valid" property value...
    goodvals = (validrange(1) <= dat & dat <= validrange(2));

    goodidx = goodhours & goodvals;

    % Build an array of daily averages for the requested field
    % (??? Skip for now - CLIPS code doesn't check values per-time either!)


    %
    % Categorize this property's values into annual and seasonal
    % ranges for "high", "very high", "somewhat low", etc., etc.
    %

    idx = goodidx;
    annual = derive_prop_ranges(dat(idx));

    % NOTE: In ORIGINAL CLIPS code, seasons started on julian days as
    % follows: spring=60, summer=120, fall=240(!), winter=330... Per
    % Jim's request this is now just based on equinoces and solstices!
    % And season names were changed to apply to N and S hemispheres.

    jasidx = (172 < jday & jday <= 264);
    idx = jasidx & goodidx;
    jas = derive_prop_ranges(dat(idx));

    ondidx = (264 < jday & jday <= 355);
    idx = ondidx & goodidx;
    ond = derive_prop_ranges(dat(idx));

    jfmidx = (355 < jday | jday <= 80);
    idx = jfmidx & goodidx;
    jfm = derive_prop_ranges(dat(idx));

    amjidx = (80 < jday & jday <= 172);
    idx = amjidx & goodidx;
    amj = derive_prop_ranges(dat(idx));

return;
