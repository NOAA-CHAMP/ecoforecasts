function stn = station_bulk_longwave(stn,afld,qfld,pfld,dsrfld,sfld,cfld,dlrf,ulrf,lrf)
%function stn = station_bulk_longwave(stn,afld,qfld,pfld,dsrfld,sfld,cfld,dlrf,ulrf,lrf)
%
% Calculate down- and upward longwave radiative fluxes using bulk formulae
% based on in situ (or other) fields in STN for air temperature STN.(AFLD),
% barometric pressure PFLD, specific humidity QFLD, insolation OR cloud cover
% DSRFLD, and sea temperature SFLD. New fields added to struct STN are total
% cloud percentage STN.(CFLD) (if necessary), and downward (DLRF), upward
% (ULRF), and net (LRF) long-wave radiative flux [W/m^2]. If new fieldnames
% are not specified, the DEFAULTS are: STN.bulk_cloud_cover, STN.bulk_dlrf,
% STN.bulk_ulrf, and STN.bulk_lrf. NOTE: If name DSRFLD == CFLD, then field
% STN.(DSRFLD) is assumed to *already contain* cloud-cover percentage.
%
% NOTE: Unlike other Ecoforecasts air-sea fluxes, ULRF is *POSITIVE UPWARD*.
%
% Example:
%  >> srvi2=load_station_data('sriv2');
%  >> srvi2=station_relhumid_to_spechumid(srvi2,'wxt_air_t','wxt_relhumid','wxt_spechumid');
%  >> srvi2=station_bulk_longwave(srvi2,'wxt_air_t','wxt_spechumid','wxt_barom', ...
%       'bic_surf_par','ctd_shallow_seatemp');
%  >> multiplot_station(srvi2,{'bulk_cloud_cover','bulk_dlrf','bulk_ulrf','bulk_lrf'});
%
% Last Saved Time-stamp: <Tue 2012-03-20 11:36:22  Lew.Gramer>

  if (~exist('cfld','var') || isempty(cfld)); cfld = 'bulk_cloud_cover'; end;
  if (~exist('dlrf','var') || isempty(dlrf)); dlrf = 'bulk_dlrf'; end;
  if (~exist('ulrf','var') || isempty(ulrf)); ulrf = 'bulk_ulrf'; end;
  if (~exist('lrf','var') || isempty(lrf)); lrf = 'bulk_lrf'; end;

  if ( isfield(stn,dlrf) ); stn = rmfield(stn,dlrf); end;
  if ( isfield(stn,ulrf) ); stn = rmfield(stn,ulrf); end;
  if ( isfield(stn,lrf) ); stn = rmfield(stn,lrf); end;

  if ( ~isfield(stn,'lon') || ~isfield(stn,'lat') )
    [stn.lon,stn.lat]=get_station_coords(stn.station_name);
  end;

  % Calculation total cloud fraction (if it was not supplied in DSRFLD)
  if ( ~strcmp(dsrfld,cfld) )
    %DEBUG:
    disp('Bulk calculating Total Cloud Cover');
    if ( isfield(stn,cfld) ); stn = rmfield(stn,cfld); end;

    [yr,mo,dy,hr] = datevec(stn.(dsrfld).date);
    yd = stn.(dsrfld).date - datenum(yr,1,1);
    [theta,cscdsr] = soradna1(yd,yr,-stn.lon,stn.lat);
    % This method produces unreliable results near sunrise/set
    eveningix = find(theta < 10);

    peakdsr = nanmax(stn.(dsrfld).data);
    Cdts = stn.(dsrfld).date;
    Craw = 100 - (100 .* stn.(dsrfld).data ./ (sind(theta).*peakdsr));

    % This method produces NO RESULT during night hours
    Cdts(eveningix) = [];
    Craw(eveningix) = [];
    Craw(Craw < 0) = 0;
    Craw(isnan(Craw)) = 0;

    % So only calculate when we can - and INTEPROLATE fluxes below...
    stn.(cfld).date = Cdts;
    stn.(cfld).data = Craw;
  end;

  dtTa = stn.(afld).date;
  [aidx,qidx] = intersect_dates(dtTa,stn.(qfld).date);  dtTa = dtTa(aidx);
  [aidx,pidx] = intersect_dates(dtTa,stn.(pfld).date);  dtTa = dtTa(aidx);
  [aidx,cidx] = intersect_dates(dtTa,stn.(cfld).date);  dtTa = dtTa(aidx);
  [aidx,sidx] = intersect_dates(dtTa,stn.(sfld).date);  dtTa = dtTa(aidx);

  [ig,ix] = intersect_dates(dtTa,stn.(qfld).date(qidx));  qidx = qidx(ix);
  [ig,ix] = intersect_dates(dtTa,stn.(pfld).date(pidx));  pidx = pidx(ix);
  [ig,ix] = intersect_dates(dtTa,stn.(cfld).date(cidx));  cidx = cidx(ix);
  [ig,ix] = intersect_dates(dtTa,stn.(sfld).date(sidx));  sidx = sidx(ix);


  stn.(lrf).date = dtTa;
  stn.(dlrf).date = dtTa;
  stn.(ulrf).date = dtTa;


  Ta = stn.(afld).data(ismember(stn.(afld).date,dtTa));
  q  = stn.(qfld).data(qidx);			% Specific Humidity [kg/kg]
  p  = stn.(pfld).data(pidx);			% Air pressure [mbar]
  if ( abs(100-nanmax(stn.(cfld).data))<50 )
    warning('Converting cloud percentage to fraction');
    C  = (stn.(cfld).data(cidx) ./ 100);	% Cloud cover [0,100]->[0.0,1.0]
  else
    C  = stn.(cfld).data(cidx);			% Cloud fraction [0.0,1.0]
  end;
  Ts = stn.(sfld).data(sidx);

  % Formulae below summarized in Weller et al (2008), and Fung et al (1984)

  % Was 0.96 - why?
  epsilon = 0.97;   % Ocean surface emissivity (see Anderson, 1952)
                    % Xue et al (1998) say: "0.97 after Anderson (1952) and Reed (1976)"
                    % Kraus and Businger suggest 0.98

  gamma = 0.622;    % "Molecular weight water / mol wt dry air" (Fung et al, 1984)
  sigma = 5.67e-8;  % Stefan-Boltzman constant
  Ta = Ta + 273.14; % convert to [K]
  Ts = Ts + 273.14; % convert to [K]
  % q = q .* 1e3;     % convert to [g/kg]  %%%% ??? No! Apparently not...
  ea = (q .* p) ./ gamma; % Near surface vapor pressure [mbar]

  % % Clark et al. (1974)
  % b = 0.62; % interpolated for latitude 25N
  % FC = 1 - b.*(C.^2); % Cloud correction factor

  % % Bunker (1976)
  % a = 0.61; % interpolated for latitude 25N
  % FC = 1 - a.*C; % Cloud correction factor

  % Hastenrath and Lamb (1978)
  b = 0.53;
  FC = 1 - b.*(C.^2); % Cloud correction factor


  % % Anderson (1952)
  % -(epsilon.*sigma.*(Ts.^4 - ((Ta.^4).*(0.74 + (O.0049.*ea)))).*FC);

  % % Berliand and Berliand (1952)
  % -((epsilon.*sigma.*(Ta.^4).*(0.39 - (O.05.*sqrt(ea))).*FC) + (4.*epsilon.*sigma.*(Ta.^3).*(Ts - Ta)));

  % % Efimova (1961)
  % -(epsilon.*sigma.*(Ta.^4).*(0.254 - (O.00495.*ea)).*FC);

  % Clark et al (1974), per Weller et al (2008)
  stn.(lrf).data = ...
      -((epsilon.*sigma.*(Ts.^4).*(0.39 - (0.05.*sqrt(ea))).*FC) + (4.*epsilon.*sigma.*(Ts.^3).*(Ts - Ta)));


  % First guess
  % Black-body radiation model of ULRF
  stn.(ulrf).data = (epsilon.*sigma.*(Ts.^4));
  % Upward should be positive, to match signs of TOMS GSIP, NCEP, and ERA
  stn.(dlrf).data = stn.(lrf).data + stn.(ulrf).data;

  % Slight tweaks
  stn.(ulrf).data = (epsilon.*sigma.*(Ts.^4)) - ((1-epsilon).*stn.(dlrf).data);
  stn.(dlrf).data = stn.(lrf).data + stn.(ulrf).data;

return;
