1;

%nix=299:699; tix=20:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % BAD
%nix=300:699; tix=20:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % JUST OK

%nix=300:699; tix=20:294; 
noLat = true;
tix=1:size(dat,1); 
while ( noLat )

  disp([min(tix),max(tix)]);

  noLon = true;
  nix=1:size(dat,2);
  while (noLon)
    d=dat(tix,nix);
    if ( ~checknans(d) )
      noLon = false;
    else
      if ( isodd(length(nix)) )
        nix(1) = [];
      else
        nix(end) = [];
      end;
      if ( isempty(nix) )
        break;
      end;
    end;
  end;

  if ( noLon )
    if ( isodd(length(tix)) )
      tix(1) = [];
    else
      tix(end) = [];
    end;
    if ( isempty(tix) )
      break;
    end;

  else
    noLat = false;
  end;

end;

if ( noLat )
  error('Found no minimal map??');
end;

