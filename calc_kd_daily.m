function stn = calc_kd_daily(stn,topfld,btmfld,zfld)
%function stn = calc_kd_daily(stn,topfld,btmfld,zfld)
%
% Add a field ['kd_1_day_' TOPFLD '_' BTMFLD] to station struct STN, with the
% Kd (exponential attenuation coefficient) for light between the two sensors
% in struct fields TOPFLD and BTMFLD. If optional arg ZFLD is a field name,
% use values in that field as inter-sensor distance between TOP and BTM; if
% ZFLD is a scalar, use that as a constant inter-sensor distance. If no ZFLD,
% assume TOP and BTM sensors are separated by exactly 1 meter of water.
%
% NOTE: Similar to CALC_KD, except ..._DAILY_MEAN values are calculated for
% both TOP and BTM, a .DEL_1_DAY_... is added to contain differences between
% those daily means, and .KD_1_DAY_... is then calculated from those fields.
%
% From the ICON/G2 definitions/documentation for the 'Kd' functions:
% // Used by Kd attenuation coefficient: assumes "surface" and "shallow" light sensors have
% //  a CONSTANT 1m of water between them. Iz = I0 * exp(-Kd*z), Kd = -(1/z)ln(Iz/I0)
% negln1m(val) = ( if ((val)<=0) then -9.0 else (- ln(val)) )
% // Used by Kd attenuation coefficient: assumes "shallow" and "deep" light sensors 2m apart
% //negln2m(val) = ( if (val<=0) then -9.0 else - ln((val)/2) ) - According to E. Stabenau - ?
% negln2m(val) = ( if ((val)<=0) then -9.0 else (- (ln(val)/2.0)) )
% 
% Last Saved Time-stamp: <Tue 2011-03-15 13:05:16  lew.gramer>

  stn = verify_variable(stn, topfld);
  stn = verify_variable(stn, btmfld);
  if ( ~exist('zfld','var') || isempty(zfld) )
    zfld = 1.0;
  elseif ( ischar(zfld) )
    stn = verify_variable(stn, zfld);
  elseif ( size(zfld) ~= [1 1] )
    error('Final arg ZFLD must be a field name or a scalar value!');
  end;

  tdts = stn.(topfld).date(stn.(topfld).data >= 0);
  tdat = stn.(topfld).data(stn.(topfld).data >= 0);
  bdts = stn.(btmfld).date(stn.(btmfld).data >= 0);
  bdat = stn.(btmfld).data(stn.(btmfld).data >= 0);

  % Tolerance for initial QA date matching is 5 minutes!
  tol = 5.0/(24.0*60.0);
% % %   [tix,bix] = intersect_dates(tdts, bdts, tol);
% % %   hrlykd = -log(bdat(bix) ./ tdat(tix));
% % %   goodix = find(isnan(hrlykd) | hrlykd > 0.10); %0/0 is ok, x/0 x>0 is not!
% % %   tdts = tdts(goodix); tdat = tdat(goodix);
% % %   bdts = bdts(goodix); bdat = bdat(goodix);

% %   [tix,bix] = intersect_dates(tdts, bdts, tol);
% %   goodix = find(tdat(tix) >= bdat(bix));
% %   tdts = tdts(tix(goodix)); tdat = tdat(tix(goodix));
% %   bdts = bdts(bix(goodix)); bdat = bdat(bix(goodix));

% %%%% DEBUG
% if (isfield(stn,[topfld '_daily_mean'])); stn = rmfield(stn, [topfld '_daily_mean']); end;
% if (isfield(stn,[btmfld '_daily_mean'])); stn = rmfield(stn, [btmfld '_daily_mean']); end;
% %%%% DEBUG

  newtfld = [topfld '_daily_mean'];
  if ( ~isfield(stn, newtfld) )
    stn.(newtfld).date = [fix(tdts(1)):fix(tdts(end))]';
    stn.(newtfld).data = repmat(nan, size(stn.(newtfld).date));
    for dix = 1:length(stn.(newtfld).date)
      tix = find(fix(tdts) == stn.(newtfld).date(dix));
      % Only pick days with reasonable amount of data
      if ( length(tix) >= 16 )
        stn.(newtfld).data(dix) = mean(tdat(tix));
      end;
    end;
  end;


  newbfld = [btmfld '_daily_mean'];
  if ( ~isfield(stn, newbfld) )
    stn.(newbfld).date = [fix(bdts(1)):fix(bdts(end))]';
    stn.(newbfld).data = repmat(nan, size(stn.(newbfld).date));
    for dix = 1:length(stn.(newbfld).date)
      bix = find(fix(bdts) == stn.(newbfld).date(dix));
      % Only pick days with reasonable amount of data
      if ( length(bix) >= 16 )
        stn.(newbfld).data(dix) = mean(bdat(bix));
      end;
    end;
  end;


  [tix,bix] = intersect_dates(stn.(newtfld).date, stn.(newbfld).date);
  tdts = stn.(newtfld).date(tix); tdat = stn.(newtfld).data(tix);
  bdts = stn.(newbfld).date(bix); bdat = stn.(newbfld).data(bix);

  ddts = tdts;
  ddat = tdat - bdat;
  kdts = tdts;
  kdat = repmat(nan, size(kdts));
  kix = find(tdat > 0 & bdat > 0);
  kdat(kix) = -log(bdat(kix) ./ tdat(kix));


  newdfld = ['del_1_day_' topfld '_' btmfld];
  newkfld = ['kd_1_day_' topfld '_' btmfld];
  % Make sure we start with a clean slate...
  if ( isfield(stn, newdfld) )
    stn = rmfield(stn, newdfld);
  end;
  if ( isfield(stn, newkfld) )
    stn = rmfield(stn, newkfld);
  end;
  stn.(newdfld).date = ddts;
  stn.(newdfld).data = ddat;
  stn.(newkfld).date = kdts;
  stn.(newkfld).data = repmat(nan, size(kdat));

  if ( isnumeric(zfld) )
    z = zfld;
    stn.(newkfld).data = kdat ./ z;
  else
    [nix,zix] = intersect_dates(stn.(newkfld).date, stn.(zfld).date);
    z = stn.(zfld).data;
    ix = find(z(zix) > 0);
    stn.(newkfld).data(nix(ix)) = kdat(nix(ix)) ./ z(zix(ix));
  end;

return;
