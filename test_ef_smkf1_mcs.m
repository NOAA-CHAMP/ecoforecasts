1;
% SCRIPT test_ef
%
% Load data for and test various ecoforecasts using ICON/matlab m-functions
%
% Last Saved Time-stamp: <Tue 2013-07-23 03:29:11 Eastern Daylight Time gramer>

set_more off;

doPrint = false;

% stanm = 'lsib4';
% stanm = 'srvi2';
% stanm = 'lppr1';
% stanm = 'fwyf1';
% stanm = 'mlrf1';
stanm = 'smkf1';
% stanm = 'sanf1';
% stanm = 'dryf1';

datapath = get_ecoforecasts_path('data/');
figspath = get_ecoforecasts_path('figs/');


% Make sure our station data is already loaded

if ( ~exist('station', 'var') || ~isfield(station,'station_name') || ~strcmpi(stanm,station.station_name) )
  station = get_station_from_station_name(stanm);
  station = load_all_ndbc_data(station);
end;


% Reload/rebuild factories
disp('Loading fact factories from CSV file...');
station = load_factories([datapath stanm '-ndbc-factories.csv'], station);
% Load any new factories or "tweaks" via a matlab call
if ( ~isempty(which(['build_factories_' stanm])) )
  disp('Building additional fact factories from M-file...');
  station = feval(['build_factories_' stanm], station);
end;


% % Verify that we now have the data we think we have
% if ( ~isfield(station, 'Name') || ~isfield(station.Name, 'data') || ...
%      ~iscell(station.Name.data) || any(~strcmpi(strtrim(station.Name.data), stanm)) )
%     error('Did we load the wrong MAT file, or is station name wrong??');
% end;



% Build our ecological forecast and evaluate it with this station's data

clear ef;

% % Rapid-Heat-Flux
% efnam = 'rhf';

% % Internal-Bore-Flux
% efnam = 'ibf';

% % Onshore-Flux
% efnam = 'of';

% % Thermal-Coral-Stress
% efnam = 'tcs';

% Mass-Coral-Stress
efnam = 'mcs';

ef = feval(['load_ef_' efnam]);


disp(sprintf('Evaluating ecoforecast "%s"...', ef.name));
% events = evaluate_ecoforecast(ef, station, datenum(2005,1,1), now);
% events = evaluate_ecoforecast(ef, station, datenum(2001,1,1), datenum(2006,12,31));
events = evaluate_ecoforecast(ef, station);



% Show the caller our results, in the Command Window and in plot figures


if ( isempty(events) )
    warning('No events produced!');

else

    % A 'period' is a sequence of two or more days, contiguous or nearly
    % so with each other, on which events of the same class occurred.
    fprintf(1, '%d event days, %d periods...\n', ...
            length(events.jday), length(events.period));

    [ uniqmos, sri ] = accumulate_monthly_sri(station, events);


    % Plot all antecedent time series, annotated by monthly cumulative S/RIs
    [vars, sens] = associate_variables(station, events);
    for vi = 1:length(vars)

        fh = figure;
        set(fh, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

        ah(1) = subplot('position', [0.10 0.28 0.88 0.64]);
        th = plot(station.(vars{vi}).date, station.(vars{vi}).data, '.');
        set(th, 'MarkerSize', 2)
        % ylim(ah(1), [min(-1,min(station.(vars{vi}).data)) max(station.(vars{vi}).data)]);
        ylim(ah(1), [min(station.(vars{vi}).data) max(station.(vars{vi}).data)]);

        % Plot S/RI for all events above the actual environmental data
        ah(2) = subplot('position', [0.10 0.05 0.88 0.20]);
        bar('v6', uniqmos, sri);
        ylabel('S/RI');

        % Sync time (X) axes, make subplots pretty together
        axes(ah(1)); datetick('x',12); set_datetick_cursor;
        axes(ah(2)); datetick('x',12); set_datetick_cursor;
        lims([1 2], 1) = xlim(ah(1)); lims([1 2], 2) = xlim(ah(2));
        minlim = min(lims(1,:)); maxlim = max(lims(2,:));
        xlim(ah(1), [minlim maxlim]); xlim(ah(2), [minlim maxlim]);
        %set(ah(1), 'XTickLabel', []);

        axes(ah(1));
        th = title(sprintf('%s vs. time (%s): Station %s, EF %s', ...
                           vars{vi}, sens{vi}, stanm, ef.name));
        set(th, 'Interpreter', 'none');
        set(gcf, 'Name', [ vars{vi} ': ' sens{vi} ]);
        if ( doPrint )
          print('-dpng', sprintf('%s-%s-%s.png', [figspath stanm], ef.name, vars{vi}));
        end;

    end;
    clear vars ah th lims minlim maxlim;

    [Y M D] = datevec(events.jday);
    yrs = unique(Y);
    nyrs = length(yrs);


    % Make a plot of seasonality in events per year
    clear hits;
    mos = 1:12;
    for mo = mos
        hits(mo) = (length(find(M == mo)) / nyrs);
    end;
    fh = figure;
    set(fh, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    plot(mos, hits);
    ylim([-1 max(hits)]);
    xlabel('Month #');
    ylabel('Events / year');
    th = title(sprintf('Event seasonality for "%s" at "%s"', ef.name, stanm));
    set(th, 'Interpreter', 'none');
    if ( doPrint )
      print('-dpng', sprintf('%s-%s-climatology.png', [figspath stanm], ef.name));
    end;


    % Make a plot of events by year
    clear hits;
    for iyr = 1:nyrs
        yr = yrs(iyr);
        hits(iyr) = length(find(Y == yr));
        fprintf(1, '%d hits in %d\n', hits(iyr), yr);
    end;
    fh = figure;
    set(fh, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    plot(yrs, hits);
    ylim([-1 max(hits)]);
    xlabel('Year');
    ylabel('Events');
    th = title(sprintf('Events per year for "%s" at "%s"', ef.name, stanm));
    set(th, 'Interpreter', 'none');
    if ( doPrint )
      print('-dpng', sprintf('%s-%s-per-year.png', [figspath stanm], ef.name));
    end;


    % Plot average length of periods, and length of longest period, by year
    clear maxes;
    [pY pM pD] = datevec(events.period);
    for iyr = 1:nyrs
        yr = yrs(iyr);
        findx = find(Y == yr);
        yrhits(iyr) = length(events.jday(findx)) / length(find(pY == yr));
        % Matlab... it's magic!
        mx = max(diff(find(diff(events.jday(findx)) > 2)));
        if ( ~isempty(mx) )
            maxes(iyr) = mx;
        else
            maxes(iyr) = yrhits(iyr);
        end;
    end;
    fh = figure;
    set(fh, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    plot(yrs, yrhits, yrs, maxes);
    ylim([-1 max(maxes)]);
    xlabel('Year');
    ylabel('Event avg/max duration (days)');
    th = title(sprintf('Event duration per year for "%s" at "%s"', ef.name, stanm));
    set(th, 'Interpreter', 'none');
    if ( doPrint )
      print('-dpng', sprintf('%s-%s-event-duration.png', [figspath stanm], ef.name));
    end;


    % Plot a histogram of period lengths over whole sample record
    fh = figure;
    set(fh, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    hist(diff(find(diff(events.jday) > 2)), max(maxes));
    th = title(sprintf('Event duration histogram for "%s" at "%s"', ef.name, stanm));
    set(th, 'Interpreter', 'none');
    if ( doPrint )
      print('-dpng', sprintf('%s-%s-duration-histo.png', [figspath stanm], ef.name));
    end;

end;

clear D;
clear findx;
clear iyr;
clear mo;
clear mos;
clear mx;
clear M;
clear nyrs;
clear pD;
clear pM;
clear pY;
clear th;
clear uniqmos;
clear vi;
clear yr;
clear yrhits;
clear Y;

set_more;
