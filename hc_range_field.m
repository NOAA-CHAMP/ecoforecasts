1;

% How far could a cold-water plug travel downhill from each point in a
% bathymetry field, without hitting a flat (beta<0.1%), reversing direction
% (angle difference > 90o), "diving" below 100 m, or hitting land?

if ( ~exist('basematfname','var') || isempty(basematfname) )
  basematfname = 'FRT_depth_and_beta_92m';
end;
if ( ~exist('dx','var') || isempty(dx) )
  dx = 92;
end;

if ( ~exist('h','var') )
  disp(['Loading bathymetry from ',[basematfname,'.mat']]);
  load([basematfname,'.mat']);
end;

%ang(jix+3:-1:jix-3,iix-3:iix+3);

checkin_count = 1000;

minrise = 0.0010*dx;

n = repmat(0,size(h));
hch = h;

DEBUG_point = [];
%DEBUG_point = [1295,3586]; % BNPON
%DEBUG_point = [1297,3574]; % BNPMI

us = [ -1 +0 +1 ;...
       -1 +0 +1 ;...
       -1 +0 +1 ];
vs = [ -1 -1 -1 ;...
       +0 +0 +0 ;...
       +1 +1 +1 ];

sz1 = size(h,1);
sz2 = size(h,2);

if ( ~isempty(DEBUG_point) )
  jix=DEBUG_point(1); iix=DEBUG_point(2);
  jixen = jix-20:jix+20;
  jixen(1>jixen|jixen>sz1) = [];
  iixen = iix-20:iix+20;
  iixen(1>iixen|iixen>sz2) = [];
  fmg;
  contourf(lon(iixen),lat(jixen),h(jixen,iixen));
  colorbar('Location','East');
end;

tic,
for jix=(sz1-1):-1:2
  if (mod(jix,checkin_count) == 0); disp([jix,sz2]); toc; tic; end;

  for iix=2:(sz2-1)
    njix = jix;
    niix = iix;

    jixen = njix-1:njix+1;
    jixen(1>jixen|jixen>sz1) = [];
    iixen = niix-1:niix+1;
    iixen(1>iixen|iixen>sz2) = [];
    hen = h(jixen,iixen);
    [mn,mnix] = min(hen(:));
    depthdif = (h(njix,niix) - mn); 
    % Give us a little more wiggle room on diagonals
    if ( mod(mnix,2) ); depthdif = depthdif / sqrt(2); end;

    while ( depthdif >= minrise && mn < 0 && mn > -80 )
      n(jix,iix) = n(jix,iix) + 1;
      hch(jix,iix) = h(njix,niix);
      u = us(mnix);
      v = vs(mnix);
      %DEBUG:      if ( jix == 1295 && iix == 3586 ); us(njix,niix) = u; vs(njix,niix) = v; end;
      if ( ~u && ~v )
        break;
      end;
      %DEBUG:      if ( jix==3638 && iix==587 ); quiver(lon(niix),lat(njix),us(mnix),vs(mnix),dx/111e3,'r'); end;
      if ( ~isempty(DEBUG_point) && all([jix,iix] == DEBUG_point) )
        quiver(lon(niix),lat(njix),us(mnix),vs(mnix),dx/111e3,'r'); 
      end;

      niix = niix+us(mnix);
      if ( 1 >= niix || niix >= sz2 )
        break;
      end;
      njix = njix+vs(mnix);
      if ( 1 >= njix || njix >= sz1 )
        break;
      end;
      %DEBUG:      if ( jix==3638 && iix==587 ); text(lon(niix),lat(njix),num2str(n(jix,iix))); keyboard; end;
      if ( ~isempty(DEBUG_point) && all([jix,iix] == DEBUG_point) )
        text(lon(niix),lat(njix),num2str(n(jix,iix)));
        keyboard;
      end;

      jixen = njix-1:njix+1;
      jixen(1>jixen|jixen>sz1) = [];
      iixen = niix-1:niix+1;
      iixen(1>iixen|iixen>sz2) = [];
      hen = h(jixen,iixen);
      [mn,mnix] = min(hen(:));
      depthdif = (h(njix,niix) - mn); 
      % Give us a little more wiggle room on diagonals
      if ( mod(mnix,2) ); depthdif = depthdif / sqrt(2); end;
      if ( mod(n(jix,iix),checkin_count) == 0 ); keyboard; end;
    end; %while ( depthdif >= minrise && mn < 0 && mn > -80 )

  end; %for jix=(sz1-1):-1:2

end; %for iix=2:(sz2-1)

hcrng = n .* dx;

toc,

if ( exist('ang','var') )
  save([basematfname,'_hc_range_depth.mat'],'lon','lat','h','bet','ang','hcrng','hch');
else
  save([basematfname,'_hc_range_depth.mat'],'lon','lat','h','bet','hcrng','hch');
end;
