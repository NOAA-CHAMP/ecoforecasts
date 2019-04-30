function stn = visedfac(stn,ef_or_facnames,resvar)
%function stn = visedfac(stn,ef_or_facnames,resvar)
%
% Visually Edit Fact-Factories for Ecoforecast
%
% Plot data fields from station struct STN, that appear as variables in any
% fact factories in STN.factories that are referenced by ecoforecast struct
% EF; red lines are plotted at the Fuzzy Lower Bound for each fuzzy in each
% factory. The user may grab any red line and move it - the resulting change
% may be stored in that factory when STN is returned. If optional arg RESVAR
% is given, it names an ecological response variable (e.g., %bleached) being
% modeled by EF: plot STN.(RESVAR) vs. each physical variable (using PLOTYY,
% qv.). If arg EF_OR_FACNAMES is a cellstr (or single string) rather than an
% ecoforecast struct, only plot variables for the named factories. DEFAULT:
% If EF_OR_FACNAMES is absent or is empty, edit all variables named in *all
% factories* in STN.factories. If any of the variables to be plotted do not
% exist in STN, use VERIFY_VARIABLE (qv.) to try to create those fields; if
% that fails, warn user but continue plotting other variables as possible. To
% finish editing, click "OK" to save factory edits, click "Cancel" or close
% figure window to discard changes in that factory, or hit CTRL-C in MATLAB
% Command Window to close all remaining figures and discard their changes.
%
% EXAMPLE #1:
%  stn = load_station_data('mlrf1');   % Load station physical data
%  ef = load_ef_mcb;                   % Load Mass-Coral-Bleaching EF
%  stn = load_frrp_bleaching(stn);     % Load TNC-FRRP bleaching data
%  stn = visedfac(stn,ef,'bleaching'); % Tweak MLRF1 factories to improve EF
%
% EXAMPLE #2:
%  stn.factories.seatemp_variability = ...
%       create_prctile_factory(stn, 'ctd_shallow_seatemp');
%  stn = visedfac(stn,'seatemp_variability');
%
% Last Saved Time-stamp: <Thu 2010-02-18 11:43:30 Eastern Standard Time Lew.Gramer>

  ef = [];
  facnames = {};

  datapath = get_ecoforecasts_path('data');

  if ( exist('ef_or_facnames','var') && ~isempty(ef_or_facnames) )
    if ( isstruct(ef_or_facnames) )
      ef = ef_or_facnames;
    elseif ( iscell(ef_or_facnames) )
      facnames = ef_or_facnames;
    elseif ( ischar(ef_or_facnames) )
      facnames = { ef_or_facnames };
    else
      error('If given, EF_OR_FACNAMES must be EF struct, or FACNAMES cellstr or name string');
    end;
  end;

  if ( ~exist('resvar','var') )
    resvar = [];
  end;


  stnm = 'Station';
  if ( isfield(stn,'station_name') )
    stnm = stn.station_name;
  elseif ( isfield(stn,'name') )
    stnm = stn.name;
  elseif ( isfield(stn,'code') )
    stnm = stn.code;
  end;

  % Reload/rebuild factories if necessary
  if ( ~isfield(stn, 'factories') )
    stn.factories = [];
    factoryfname = fullfile(datapath, [stnm '-factories.csv']);
    if ( exist(factoryfname,'file') )
      fprintf(1, 'Loading fact factories from CSV file:\n "%s"\n', factoryfname);
      stn = load_factories(factoryfname, stn);
    end;
    % Load any new factories or "tweaks" via a matlab call
    factorymfile = which(['build_factories_' stnm]);
    if ( ~isempty(factorymfile) )
      fprintf(1, 'Building additional fact factories from M-file:\n "%s"\n', factorymfile);
      stn = feval(['build_factories_' stnm], stn);
    end;
  end;


  if ( ~isempty(ef) )
    rules = fieldnames(ef.rules);
    for rulix = 1:length(rules)
      facnames = { facnames{:} ef.rules.(rules{rulix}).facts{:} };
    end;
    facnames = unique(facnames);
  else
    if ( isempty(facnames) )
      if ( isfield(stn, 'factories') && isstruct(stn.factories) )
        facnames = fieldnames(stn.factories);
      end;
    end;
  end;

  % User didn't provide factory names and we can't make any, so give up
  if ( isempty(facnames) )
    error('No factory-names could be found in args, EF, or STN.factories!');
  end;


  % If STN has no factories, create ugly new ones for the user to edit
  if ( ~isfield(stn, 'factories') || ~isstruct(stn.factories) )
    badix = [];
    for facix = 1:length(facnames)
      facname = facnames{facix};
      varname = default_factory_var(facname);
      if ( isempty(varname) )
        warning('Skipping factory "%s": No default variable name!', facname);
        badix(end+1) = facix;
      else
        stn = verify_variable(stn, varname);
        if ( ~isfield(stn,varname) )
          warning('Skipping factory "%s": Cannot create default variable %s', facname, varname);
          badix(end+1) = facix;
        else
          % This call creates a very rough set of fuzzy ranges for the new
          % factory, based on simple percentiles of the data distribution
          stn.factories.(facname) = create_prctile_factory(stn,varname);
        end;
      end;
    end;
    facnames(badix) = [];
  end;


  %%%% EARLY RETURN
  %%%% EARLY RETURN
  if ( isempty(facnames) )
    warning('None of requested factory-names can be loaded or created! Returning early...');
    return;
  end;


  varnames = {};
  varfacnames = {};
  for facix = 1:length(facnames)
    fac = facnames{facix};
    if ( iscell(stn.factories.(fac).variable) )
      varnames = { varnames{:} stn.factories.(fac).variable{:} };
    else
      varnames = { varnames{:} stn.factories.(fac).variable };
    end;
    varfacnames = { varfacnames{:} fac };
  end;
  [varnames, facix] = unique(varnames);
  varfacnames = varfacnames(facix);


  badix = [];
  for vix = 1:numel(varnames)
    stn = verify_variable(stn, varnames{vix});
    if ( ~isfield(stn,varnames{vix}) )
      warning('Unable to create missing field %s', varnames{vix});
      badix(end+1) = vix;
    end;
  end;

  varnames(badix) = [];
  varfacnames(badix) = [];
  nvars = numel(varnames);


  %%%% EARLY RETURN
  %%%% EARLY RETURN
  if ( isempty(varnames) )
    warning('No variables can be found or built for factory-names! Returning early...');
    return;
  end;



  innerlimits = [-Inf +Inf];
  outerlimits = [+Inf -Inf];
  if ( ~isempty(resvar) )
    innerlimits = stn.(resvar).date([1 end]);
    outerlimits = stn.(resvar).date([1 end]);
  end;
  for vix = 1:nvars
    innerlimits = [ max(innerlimits(1), stn.(varnames{vix}).date(1)) , ...
                    min(innerlimits(2), stn.(varnames{vix}).date(end)) ];
    % Do not attempt to adjust your television set
    outerlimits = [ min(outerlimits(1), stn.(varnames{vix}).date(1)) , ...
                    max(outerlimits(2), stn.(varnames{vix}).date(end)) ];
  end;


  resfig = [];
  % Plot response variable as its own graph, just for convenience
  if ( ~isempty(resvar) )
    resfig = figure;
    plot(stn.(resvar).date, stn.(resvar).data, 'b-d');
    xlim(gca, innerlimits);
    maxigraph;
    datetick3;
    ttl = sprintf('%s %s (for comparison)', stnm, upper(resvar));
    titlename(strrep(ttl,'_','\_'));
  end;


  for vix = 1:nvars

    var = varnames{vix};
    fac = varfacnames{vix};

    fhs(vix) = figure;

    if ( ~isempty(resvar) )
      [ix1,ix2] = intersect_dates(stn.(resvar).date, stn.(var).date, 1);
      if ( isempty(ix1) )
        warning('No %s dates match variable %s!',resvar,varnames{ix});
      end;
    end;

    if ( isempty(resvar) )
      hvar = plot( stn.(var).date, stn.(var).data );
      setappdata(gcf, 'factoryaxes', gca);

    else
      [axs,hres,hvar] = ...
          plotyy( stn.(resvar).date, stn.(resvar).data, ...
                  stn.(var).date, stn.(var).data );
      ylabel(axs(1), upper(strrep(resvar,'_','\_')));
      set(hres, 'Marker','diamond');
      set(axs(1), 'YGrid','off', 'YMinorGrid','off', 'YMinorTick','off');

      % Link X-axes in order to synchronize pan/zoom
      linkaxes(axs,'x');
      xlim(axs(1), innerlimits);

      setappdata(gcf, 'factoryaxes', axs(2));
      axes(axs(2));

    end;

    setappdata(gcf, 'limits', outerlimits);
    setappdata(gcf, 'factory', fac);

    ylabel(upper(strrep(var,'_','\_')));
    xlim(innerlimits);
    % ylim([nanmin(stn.(var).data), nanmax(stn.(var).data)]);
    set(gca, 'YTickMode','auto', 'YTickLabelMode','auto');
    set(gca, 'YGrid','on', 'YMinorGrid','on', 'YMinorTick','on');

    % Allow the "bar" for each fuzzy lower-bound to be raised or lowered
    f = stn.factories.(fac);
    for idx = 1:length(f.fuzzies)
      rng = f.fuzzies{idx};
      [hl,ht] = annotline([], rng{2}(1), ['\uparrow' rng(1)], 'red');
      setappdata(hl,'htext',ht);
      set(hl,'ButtonDownFcn',@visedfac_fuzzy_btndn);
      hlines(idx) = hl;
    end;
    [hl,ht] = annotline([], rng{2}(2), ['\uparrow' 'unbelievably-high'], 'red');

    ylim([f.fuzzies{1}{2}(1) f.fuzzies{end}{2}(2)]);

    setappdata(hl,'htext',ht);
    set(hl,'ButtonDownFcn',@visedfac_fuzzy_btndn);
    hlines(idx+1) = hl;

    setappdata(gcf,'hlines',hlines);

    maximize_graph;
    datetick3;

    ttl = sprintf('%s %s (Factory %s)', stnm, upper(var), fac);
    if ( ~isempty(resvar) )
      ttl = [ ttl ' vs. ' upper(resvar) ];
    end;
    titlename(strrep(ttl,'_','\_'));

    okic = uicontrol('Style', 'pushbutton', 'String', 'Confirm', ...
                     'Units', 'normalized', 'Value', 1, 'Parent', gcf, ...
                     'Pos', [0.01 0.81 0.06 0.08]);
    set(okic, 'Callback', {@visedfac_confirm_cancel, gcf, true});
    cxlic = uicontrol('Style', 'pushbutton', 'String', 'Cancel', ...
                      'Units', 'normalized', 'Value', 1, 'Parent', gcf, ...
                      'Pos', [0.01 0.71 0.06 0.08]);
    set(cxlic, 'Callback', {@visedfac_confirm_cancel, gcf, false});

    % Make sure Figure toolbar is showing...
    set(gcf,'Toolbar','figure');

  end; % for vix

  try
    for fi = 1:length(fhs)
      % Block on EACH of the figures in turn: if the user closes all of them,
      % or hits CTRL-C (which will error out of this loop), only then return.
      if ( ishandle(fhs(fi)) )
        %%%% 
        %%%% UIWAIT - blocking call
        %%%% 
        if ( strcmp(get(fhs(fi),'Visible'),'on') )
          uiwait(fhs(fi));
        end;

        OK = [];
        if ( ishandle(fhs(fi)) )
          OK = getappdata(fhs(fi), 'confirmed');
        end;

        if ( ~isempty(OK) && OK && isappdata(fhs(fi), 'factory') )
          fac = getappdata(fhs(fi), 'factory');
          hlines = getappdata(fhs(fi), 'hlines');

          f = stn.factories.(fac);
          for idx = 1:length(f.fuzzies)
            rng = f.fuzzies{idx};

            ys = get(hlines(idx),'YData');
            rng{2}(1) = ys(1);

            ys = get(hlines(idx+1),'YData');
            rng{2}(2) = ys(1);

            stn.factories.(fac).fuzzies{idx} = rng;
          end;
          %DEBUG:          fprintf(1, 'Saved Changes: "stn.factories.%s"\n', fac);

        else
          if ( ishandle(fhs(fi)) )
            fhname = get(fhs(fi),'Name');
          else
            fhname = sprintf('Figure %g', fhs(fi));
          end;
          fprintf(1, 'CANCELLED: "%s"\n', fhname);
        end;

      end;
    end;

  catch
    warning('CAUGHT INTERRUPT: Some factory changes may have been lost!');
  end;

  for fi = 1:length(fhs)
    if ( ishandle(fhs(fi)) )
      close(fhs(fi));
    end;
  end;

  if ( ishandle(resfig) )
    close(resfig);
  end;

return;



%%%%%%%%%% 
%%%%%%%%%% PRIVATE FUNCTIONS
%%%%%%%%%% 

function visedfac_confirm_cancel(uich, eventdata, fh, confirmed)
%function visedfac_confirm_cancel(uich, eventdata, fh, confirmed)
% Close figure window, either saving changes, or canceling

  if ( confirmed )
    setappdata(fh, 'confirmed', 1);
  end;
  %%%% 
  %%%% UIRESUME - unblocks uiwait call in primary function above
  %%%% 
  uiresume(fh);

  % Either way, user is done looking at this figure now...
  if ( ishandle(fh) )
    set(fh, 'Visible', 'off');
  end;

return;


function visedfac_fuzzy_btndn(src, eventdata)
% Let the user move a fuzzy lower-bound on the graph

  fh = gcbf;
  ah = getappdata(fh, 'factoryaxes');
  fac = getappdata(fh, 'factory');
  limits = getappdata(fh, 'limits');

  hl = gcbo;
  ht = getappdata(hl,'htext');
  fuzzy = get(ht,'string');

  %DEBUG:  ['Editing ' fac '->' fuzzy],

  point1 = get(ah, 'CurrentPoint');    % button down detected
  begval = point1(1,2);
  %DEBUG:  begval,

  finalRect = rbbox;      % return (idiotically) plotted units

  point1 = get(ah, 'CurrentPoint');    % button up
  endval = point1(1,2);
  %DEBUG:  endval,

  set(hl,'YData',[endval endval]);
  tpos = get(ht,'Position');
  set(ht,'Position',[tpos(1) endval tpos(3)]);

return;
