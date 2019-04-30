 1;

 more_status = get(0, 'More'); more('off');

 % Retrieve and store data in a path relative to this M-file's local directory
 if ( ~exist('datapath', 'var') || isempty(datapath) )
   datapath = get_ecoforecasts_path('data');
 end;
 if ( ~exist('figspath', 'var') || isempty(figspath) )
   figspath = get_ecoforecasts_path('figs');
 end;

 period = 'annual';
 cutoff_prctile = 93;

 onehour = (1.0/24.0);
 mins29 = (29.0/(24.0*60.0));

 % stnam = 'lkwf1';
 % stnam = 'fwyf1';
 % stnam = 'mlrf1';
 % stnam = 'lonf1';
 % stnam = 'smkf1';
 % stnam = 'sanf1';
 % stnam = 'dryf1';
 stnam = '42003';

 if ( ~exist('station', 'var') || isempty(station) )
   station = load_all_ndbc_data([], stnam);
 end;

 % Sanity check our station data
 if (~isfield(station,'station_name') || ~strcmp(station.station_name,stnam))
   error('Did we load the wrong data?! Station_name ~= "%s"', stnam);
 end;


 % Get ready to save our intermediate results as we go along...
 statemfname = fullfile(datapath, ...
                        sprintf('SOM-%s-%s-%g.mat', ...
                                stnam, period, cutoff_prctile));

 if ( ~exist('ndays', 'var') || isempty(ndays) )
   ndays = 14; % Days per "frame"
 end;
 daysbefore = ceil(ndays / 2);
 nhrs = (ndays*24);

 if ( ~exist('mapdims', 'var') || isempty(mapdims) )
   mapdims = [3 3];
 end;
 nmaps = mapdims(1) * mapdims(2);

 anfld = 'sea_t_1_day_deviation_3_day_average';


 ix = 0;

 % Analyze air-sea flux variables
 % ix=ix+1; fcfld{ix} = 'wind_stress';
 % ix=ix+1; fcfld{ix} = 'wind_stress_1_day_integral';
 % ix=ix+1; fcfld{ix} = 'wind_stress_u_7_day_deviation_sum_wind_stress_v';
 % ix=ix+1; fcfld{ix} = 'net_heat_flux';
 % ix=ix+1; fcfld{ix} = 'net_heat_flux_1_day_integral';

 % *OR* analyze atmospheric variability metrics
 ix=ix+1; fcfld{ix} = 'wind1_speed_3_day_average';
 ix=ix+1; fcfld{ix} = 'wind1_u_7_day_deviation_sum_wind1_v';
 ix=ix+1; fcfld{ix} = 'air_t_1_day_deviation_3_day_average';
 ix=ix+1; fcfld{ix} = 'barom_1_day_deviation_3_day_average';

 nvars = length(fcfld);

 station = verify_variable(station, anfld);
 for ix = 1:nvars
   station = verify_variable(station, fcfld{ix});
 end;

 % Figure smarter way to do this later on...
 station = blank_derived_var(station, 'sea_t', anfld, [], 3);
 station = blank_derived_var(station, 'wind1_speed', fcfld{1}, [], 3);
 station = blank_derived_var(station, 'wind1_speed', fcfld{2}, [], 7);
 station = blank_derived_var(station, 'air_t', fcfld{3}, [], 3);
 station = blank_derived_var(station, 'barom', fcfld{4}, [], 3);

 % LATER... add code to LOAD from our MAT file if we already saved one!

 fprintf('Building %d-day frames around %g%% events...\n', ndays, cutoff_prctile);

 [anix,fcix1] = intersect_dates(station.(anfld).date, station.(fcfld{1}).date);

 andts = station.(anfld).date(anix);
 andat = station.(anfld).data(anix);

 anmn = nanmean(andat);
 ansd = nanstd(andat);


 dv = datevec(andts([1 end]));
 begyr = dv(1,1);
 endyr = dv(2,1);

 ix = 1;
 fcdts{1} = station.(fcfld{ix}).date(fcix1);
 fcdat{1} = station.(fcfld{ix}).data(fcix1);
 fcmn(1) = nanmean(station.(fcfld{ix}).data(fcix1));
 fcsd(1) = nanstd(station.(fcfld{ix}).data(fcix1));
 for ix = 2:nvars
   [ig,fcix2] = intersect_dates(station.(fcfld{1}).date(fcix1), station.(fcfld{ix}).date);
   fcdts{ix} = station.(fcfld{ix}).date(fcix2);
   fcdat{ix} = station.(fcfld{ix}).data(fcix2);
   fcmn(ix) = nanmean(station.(fcfld{ix}).data(fcix2));
   fcsd(ix) = nanstd(station.(fcfld{ix}).data(fcix2));
 end;


%  % Normalize all our variables
%  andat = (andat - anmn) ./ ansd;
%  for ix = 1:nvars
%    fcdat{ix} = (fcdat{ix} - fcmn(ix)) ./ fcsd(ix);
%  end;


 % Find all 'events' - periods when response variable is above cutoff value
 cutoff = prctile(andat, cutoff_prctile);

 ix = find(andat > cutoff);

 dd = diff(andts(ix));
 evtixix = find(dd > (ndays+1));
 evtix = ix(evtixix);
 evtdts = andts(evtix);

 interpmthd = 'linear';
 extrapval = nan;

 dts = andts(1):onehour:andts(end);
 andat = interp1(andts, andat, dts, interpmthd, extrapval);
 andts = dts;
 for ix = 1:nvars
   fcdat{ix} = interp1(fcdts{ix}, fcdat{ix}, dts, interpmthd, extrapval);
   fcdts{ix} = dts;
 end;

 allfdts = [];
 allfdat = [];
 for ix = 1:length(evtdts)

   baddata = false;

   % Build a data row-vector by concatenating the two-week interpolated
   % samples from the response variable and each of the forcing variables

   evtix = find( (abs(andts - evtdts(ix)) <= mins29), 1 );
   if ( 0 > evtix-(daysbefore*24) || evtix+(daysbefore*24)-1 > numel(andts) )
     warning('Skipping event %d - too near limits of %s??', ix, anfld);
     continue;
   end;
   fdt = andts((evtix-(daysbefore*24)):(evtix+(daysbefore*24)-1));
   fda = andat((evtix-(daysbefore*24)):(evtix+(daysbefore*24)-1));
   if ( any(isnan(fda)) )
     baddata = true;
     warning('Skipping event %d - bad data in %s??', ix, anfld);
   end;
   fdts = fdt;
   fdat = fda;

   for fcix = 1:nvars
     evtix = find( (abs(fcdts{fcix} - evtdts(ix)) <= mins29), 1 );
     fdt = fcdts{fcix}((evtix-(daysbefore*24)):(evtix+(daysbefore*24)-1));
     fda = fcdat{fcix}((evtix-(daysbefore*24)):(evtix+(daysbefore*24)-1));
     if ( ~baddata && any(isnan(fda)) )
       baddata = true;
       warning('Skipping event %d - bad data in %s??', ix, fcfld{fcix});
     end;
     fdts = [ fdts , fdt ];
     fdat = [ fdat , fda ];
   end;

   % Add current "event" row-vector to the training data matrix
   % (But only if we have clean data for all variables!)
   if ( ~baddata )
     allfdts = [ allfdts ; fdts ];
     allfdat = [ allfdat ; fdat ];
   end;

 end;


 % Once we reach this point - save our state...
 fprintf('Saving state...\n');
 save(statemfname);


 % LATER... add code to LOAD from our MAT file if we already created one!


 %
 % Self-Organizing Map (Artificial Neural Network) training and analysis
 %

 fprintf('Training Self-Organizing Map...\n');


 sominit = 'randinit';
 somneigh = 'ep';
 somneigh = 'gaussian';
 sm = som_make(allfdat, 'msize',mapdims, sominit, 'batch', 'tracking',0, ...
               'shape','sheet', 'training','long', 'neigh',somneigh);

 % Store various other run data in the 'sm' struct
 sm.station_name = stnam;
 sm.period = period;
 sm.cutoff_prctile = cutoff_prctile;
 sm.cutoff = cutoff;

 % Add Best-Matching Units vector as a member of each 'sm' struct.
 % NOTE: For SOM Unit (map-node, "mode") 'n', framedts((sm.bmus==n),[1 end])
 % are the start- and end-dates for all the frames best matched by that Unit.
 [sm.bmus, sm.qerrs] = som_bmus(sm, allfdat);

 % Calculate number of matched sample frames for each "mode" (Unit).
 sm.permatched = som_hits(sm, allfdat);
 sm.pctmatched = sm.permatched * (100 / size(allfdat,1));

 % Sort BMUs by number of frames, from most-matched to least.
 [ign, sm.mode_order] = sort(sm.permatched, 'descend');


%  % De-normalize the variables, and each unit ("mode") in the SOM
%  andat = (andat .* ansd) + anmn;
%  for ix = 1:nvars
%    fcdat{ix} = (fcdat{ix} .* fcsd(ix)) + fcmn(ix);
%  end;
%  for cbix = 1:nmaps
%    cb = reshape(sm.codebook(cbix,:), [nhrs (nvars+1)]);
%    cb(:,1) = (cb(:,1) .* ansd) + anmn;
%    for ix = 1:nvars
%      cb(:,ix+1) = (cb(:,ix+1) .* fcsd(ix)) + fcmn(ix);
%    end;
%    sm.codebook(cbix,:) = cb(:)';
%  end;


 % Once we reach this point - save our state AGAIN...
 fprintf('Saving state...\n');
 save(statemfname);


 figure;

 set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
 set(gcf, 'DefaultAxesFontSize', [8]);
 cmins = min(sm.codebook); cmn = min(reshape(cmins, [nhrs (nvars+1)]));
 cmaxs = max(sm.codebook); cmx = max(reshape(cmaxs, [nhrs (nvars+1)]));
 cmin = min([-abs(cmn) ; -abs(cmx)]); cmax = max([abs(cmn) ; abs(cmx)]);
 cmin = cmn; cmax = cmx;
 %%%% ???
 cmin = repmat(0, size(cmn)); cmax = cmx;
 for orderix = 1:length(sm.mode_order)
   % Plot "modes" from most frequently matched to least
   ix = sm.mode_order(orderix);
   subplot(mapdims(2), mapdims(1), orderix);
   hold on;
   x = 1:nhrs;
   cdbk = reshape(sm.codebook(ix,:), [nhrs (nvars+1)]);
   do_tvar_peak_plots(x, cdbk, cmin, cmax);
   hold off;
   title(sprintf('#%d (%d frames - %0.1f%%)', ...
                 ix, sm.permatched(ix), sm.pctmatched(ix)));
 end;
 pltclrs = 'bgrcmykw';
 fcfldlbl = '';
 for ix = 1:nvars
   fcfldlbl = [ fcfldlbl strrep(lower(fcfld{ix}),'_','\_')  ' (' pltclrs(ix+1) '), '];
 end;
 suplabel(fcfldlbl, 'x');
 suplabel(sprintf( 'SOM "Modes" %s %d-%d %s Events (>%g%%,%g), %dH frames (N=%d): %s (%s)', ...
                   upper(stnam), begyr, endyr, upper(period), cutoff_prctile, cutoff, ...
                   nhrs, size(allfdat,1), strrep(lower(anfld),'_','\_'), pltclrs(1) ), 't');
 pfname = sprintf('SOM-%s-%s-%s%g-%02d-%02d', lower(stnam), ...
                  period, 'pct', cutoff_prctile, mapdims(1), mapdims(2));
 print('-dtiff', '-r300', fullfile(figspath, [pfname '.tiff']));

 % Done
 more(more_status);
