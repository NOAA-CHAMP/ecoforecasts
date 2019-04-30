function stn = calc_kd(stn, topfld, btmfld, zfld)
%function stn = calc_kd(stn, topfld, btmfld, zfld)
%
% Add a field ['kd_' TOPFLD '_' BTMFLD] to station struct STN, containing the
% Kd (exponential attenuation coefficient) for light between the two sensors
% in struct fields TOPFLD and BTMFLD. If optional arg ZFLD is a field name,
% use values in that field as inter-sensor distance between TOP and BTM; if
% ZFLD is a scalar, use that value as a constant inter-sensor distance. If
% ZFLD is not given, assume sensors TOP and BTM are separated by 1 meter.
%
% From the ICON/G2 definitions/documentation for the 'Kd' functions:
% // Used by Kd attenuation coefficient: assumes "surface" and "shallow" light sensors have
% //  a CONSTANT 1m of water between them. Iz = I0 * exp(-Kd*z), Kd = -(1/z)ln(Iz/I0)
% negln1m(val) = ( if ((val)<=0) then -9.0 else (- ln(val)) )
% // Used by Kd attenuation coefficient: assumes "shallow" and "deep" light sensors 2m apart
% //negln2m(val) = ( if (val<=0) then -9.0 else - ln((val)/2) ) - According to E. Stabenau - ?
% negln2m(val) = ( if ((val)<=0) then -9.0 else (- (ln(val)/2.0)) )
% 
% Last Saved Time-stamp: <Wed 2010-06-16 13:47:17 Eastern Daylight Time Lew.Gramer>

  stn = verify_variable(stn, topfld);
  stn = verify_variable(stn, btmfld);
  if ( ~exist('zfld','var') || isempty(zfld) )
    zfld = 1.0;
  elseif ( ischar(zfld) )
    stn = verify_variable(stn, zfld);
  elseif ( size(zfld) ~= [1 1] )
    error('Final arg ZFLD must be a field name or a scalar value!');
  end;

  [tix,bix] = intersect_dates(stn.(topfld).date, stn.(btmfld).date);
  tdat = stn.(topfld).data(tix);
  bdat = stn.(btmfld).data(bix);

  tdat(tdat <= 0) = nan;
  bdat(bdat <= 0) = nan;

  newfld = ['kd_' topfld '_' btmfld];
  % Make sure we start with a clean slate...
  if ( isfield(stn, newfld) )
    stn = rmfield(stn, newfld);
  end;
  stn.(newfld).date = stn.(topfld).date(tix);
  stn.(newfld).data = repmat(-9.0, size(stn.(newfld).date));

  ix = find(tdat ~= 0);
  val = stn.(newfld).data;
  val(ix) = bdat(ix) ./ tdat(ix);

  if ( isnumeric(zfld) )
    disp(['Using scalar intersensor distance ' num2str(zfld) 'm']);
    nix = 1:length(stn.(newfld).data);
    z = repmat(zfld, size(stn.(newfld).data));
  else
    [nix,zix] = intersect_dates(stn.(newfld).date, stn.(zfld).date);
    % z = stn.(zfld).data(zix)';
    z = stn.(zfld).data(zix);
  end;

  % Only use physically possible values!
  ix = find(val(nix) > 0);

  sun_angle_correction = 1.0;
  if ( isfield(stn,'lon') )
    [yds,yrs] = get_yearday(stn.(newfld).date);
    [theta,ig] = soradna1(yds,yrs,-stn.lon,stn.lat);
    sun_angle_correction = secd(90-theta(ix));
  end;

  stn.(newfld).data(nix(ix)) = -log(val(nix(ix))) ./ (z(ix) .* sun_angle_correction);

  stn.(newfld).data(-0.1 > stn.(newfld).data | stn.(newfld).data > 5) = nan;

return;
