1;

spp_pairs=[]; mind=[]; minpix=[]; minpst={};

for ix = 1:numel(res.ms)
  minmst{ix} = res.ms(ix).station_name;
  minpst{ix} = '';
  mind(ix) = +inf;
  minpix(ix) = 0;

  for stix = 1:numel(res.ps)
    d = station_dist(res.ms(ix),res.ps(stix));
    if ( d < mind(ix) )
      mind(ix) = d;
      minpix(ix) = stix;
      minpst{ix} = res.ps(stix).station_name;
    end;
  end;

  if ( mind(ix) > eps )
    if ( mind(ix) < 1 )
      warning('M. cav. site %s: Nearest P. ast. site %s %g km away', ...
              res.ms(ix).station_name,res.ps(minpix(ix)).station_name,mind(ix));
    else
      error('No P. ast. site anywhere near M. cav. site %s',res.ms(ix).station_name);
    end;
  end;

  spp_pairs(ix).mix = ix;
  spp_pairs(ix).pix = minpix(ix);
  spp_pairs(ix).mnm = minmst{ix};
  spp_pairs(ix).pnm = minpst{ix};
  spp_pairs(ix).m = res.ms(ix);
  spp_pairs(ix).p = res.ps(minpix(ix));
end;

fprintf(1,'Species site mismatches:\nP. ast.\tM. cav.\t Dist [m] \n==========================\n');
for ix = 1:numel(res.ms)
  if ( mind(ix) > eps || ~strcmpi(minpst{ix},res.ms(ix).station_name) )
    fprintf(1,'%s \t %s \t %g \n',minpst{ix},res.ms(ix).station_name,round(mind(ix)*1e3));
  end;
end;

clear ix stix d
