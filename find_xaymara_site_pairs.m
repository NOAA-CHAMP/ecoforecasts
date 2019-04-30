1;

for ix = 1:numel(res.ps)
  x0 = res.ps(ix);
  for stix = 1:numel(res.ps)
    if ( stix > ix )
      find_xaymara_site_pairs_filter(res.ps(ix),res.ps(stix),'P')
    end;
  end;
end;

fprintf(1,'\n');

for ix = 1:numel(res.ms)
  for stix = 1:numel(res.ms)
    if ( stix > ix )
      find_xaymara_site_pairs_filter(res.ms(ix),res.ms(stix),'M')
    end;
  end;
end;

clear ans ax d ix stix x x0


%%%% INTERNAL FUNCTION

function find_xaymara_site_pairs_filter(x0,x,typ,pct,mind,mindbeta)
%function find_xaymara_site_pairs_filter(x0,x,typ,pct,mind,mindbeta)
% Print a one-line report if site X (STRUCT) is within MIND km (DEFAULT: 10)
% of site X0, has depth within PCT of X0's (DEFAULT: 0.25), and has seafloor
% slope at least MINDBETA (DEFAULT: 0.01) different from X0's beta.

  if ( ~exist('pct','var') || isempty(pct) )
    pct = 0.25;
  end;
  if ( ~exist('mind','var') || isempty(mind) )
    mind = 10.0;
  end;
  if ( ~exist('mindbeta','var') || isempty(mindbeta) )
    mindbeta = 0.010;
  end;
  d = station_dist(x, x0);
  dnd = abs(x.ngdc_depth - x0.ngdc_depth);
  dvd = abs(x.very_hires_smooth_depth - x0.very_hires_smooth_depth);
  dxd = abs(x.depth - x0.depth);
  db = abs(x.very_hires_smooth_beta - x0.very_hires_smooth_beta);
  if ( d < mind && db > mindbeta && ...
       (dnd < abs(x0.ngdc_depth*pct) || dvd < abs(x0.very_hires_smooth_depth*pct) || dxd < abs(x0.depth*pct)) )
    fprintf(1,'%s: % 5s,% 5.1f,% 5.1f,% 5.1f,%.3f,%.3f  % 5s,% 5.1f,% 5.1f,% 5.1f,%.3f,%.3f  %4.0fm % 5.1f/% 5.1f %.3f % 3.0f % 3.0f\n',...
            typ,...
            x0.station_name,x0.depth,x0.ngdc_depth,x0.very_hires_smooth_depth,x0.ngdc_beta,x0.very_hires_smooth_beta,...
            x.station_name, x.depth,x.ngdc_depth,x.very_hires_smooth_depth,x.ngdc_beta,x.very_hires_smooth_beta,...
            d*1e3,dnd,dvd,db,x0.Ng,x.Ng);
  end;

end
