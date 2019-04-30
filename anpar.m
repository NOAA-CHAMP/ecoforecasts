1;

tfld = 'bic_surf_par';
bfld = 'bic_shallow_par';
dfld = 'bic_dpar1m';

if ( ~exist('srvi2','var') || isempty(srvi2) )
  srvi2 = load_station_data('srvi2');
end;

% Tolerance for initial QA date matching is 5 minutes!
tol = 5.0/(24.0*60.0);
tol = [];
[tix,bix] = intersect_dates(srvi2.(tfld).date,srvi2.(bfld).date,tol);

if (isfield(srvi2,dfld)); srvi2=rmfield(srvi2,dfld); end;
srvi2.(dfld).date = srvi2.(tfld).date(tix);
srvi2.(dfld).data = srvi2.(tfld).data(tix) - srvi2.(bfld).data(bix);

multiplot_station(srvi2,{tfld,bfld,dfld},'SRVI2 \DeltaPAR_1_m QA');
