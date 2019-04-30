function stn = station_wind_to_wave(stn,wspdfld,wdirfld,wvperfld,wvhgtfld,wvdirfld,wvufld,wvvfld)
%function stn = station_wind_to_wave(stn,wspdfld,wdirfld,wvperfld,wvhgtfld,wvdirfld,wvufld,wvvfld)
%
% Add fields with estimates of wind-wave height, STN.(WVHGTFLD), and wave
% period, STN.(WVPERFLD), from station wind speed time series STN.(WSPDFLD).
% Formula from Fairall MATLAB script COR3_0AH.M. (Tp from Pierson-Moskowitz
% relation, significant wave height from one of several relations reviewed in
% Huang et al., 1990; both referenced in e.g., Bourassa et al, 2001. Taylor
% and Yelland, 2001, another formulation for Hs was also considered.) Copies
% wind direction STN.(WDIRFLD) to new wave direction STN.(WVDIRFLD) field.
% If optional WVUFLD and WVVFLD are specified, also calculates components
% STN.(WVUFLD) and STN.(WVVFLD) of peak wavenumber vector time series.
%
% NOTE: Assumes STN.(WSPDFLD) wind speed is in [kts].
%
% Last Saved Time-stamp: <Thu 2012-01-26 16:17:28  Lew.Gramer>

  [wvper,wvhgt] = wind_to_wave(stn.(wspdfld).data);

  stn.(wvperfld).date = stn.(wspdfld).date;
  stn.(wvperfld).data = wvper;

  stn.(wvhgtfld).date = stn.(wspdfld).date;
  stn.(wvhgtfld).data = wvhgt;

  stn.(wvdirfld).date = stn.(wdirfld).date;
  stn.(wvdirfld).data = stn.(wdirfld).data;

  % If caller requests it, also calculate peak wave vector components "U" and "V"
  if ( exist('wvufld','var') && exist('wvvfld','var') && ~isempty(wvufld) && ~isempty(wvvfld) )
    [Tix,Dix] = intersect_dates(stn.(wvperfld).date,stn.(wvdirfld).date);
    T = stn.(wvperfld).data(Tix);
    D = stn.(wvdirfld).data(Dix);
    % Calculate "u" and "v" peak wavenumber vector components by weighting
    % PEAKWAVEDIR with putative (deep-water) wave celerity: thus. the lower
    % frequency waves dominate when averaged in with shorter wavelengths.
    % (As for any vector, interpolating PEAKWAVEDIR means special handling!)
    if ( isfield(stn,'lat') ); lat=stn.lat; else lat=25; end;
    g = sw_g(lat,0);
    c = g .* T ./ (2*pi);
    % % *OR* just do a simple vector averaging with all weights == 1??
    % c = 1;
    [u,v] = spddir_to_uv(c,D);
    stn.(wvufld).date = stn.(wvdirfld).date(Dix);
    stn.(wvufld).data = u;
    stn.(wvvfld).date = stn.(wvdirfld).date(Dix);
    stn.(wvvfld).data = v;
  end;


return;
