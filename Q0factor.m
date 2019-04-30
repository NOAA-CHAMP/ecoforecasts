function Q0_factor = Q0factor(t,s,d)
%function Q0_factor = Q0factor(t,s,d)
%
% Calculate divisor for heat flux(es) in heat budget equation
%

  representative_mean_salinity = 36;
  if ( ~exist('s','var') || isempty(s) )
    s = representative_mean_salinity;
  end;

  representative_depth_meters = 2;
  if ( ~exist('d','var') || isempty(d) )
    d = representative_depth_meters;
  end;

  if ( numel(t) > 1 )
    t = nanmean(t(:));
  end;
  if ( numel(s) > 1 )
    s = nanmean(s(:));
  end;
  if ( numel(d) > 1 )
    d = nanmean(d(:));
  end;

  % Cp and rho from [Fofonoff and Millard, 1983]
  rho = sw_dens0( s, t );
  Cp = sw_cp( s, t, d );
  h = d;

  Q0_factor = (rho*Cp*h);

return;
