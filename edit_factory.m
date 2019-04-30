function stn = edit_factory(stn, fac)
%function stn = edit_factory(stn, fac)
%
% Pop up a series of dialog boxes that allow the user to edit a fact-factory
% named "fac" in the station struct "stn". If stn.factories or the subfield
% stn.factories.(fac) do not exist, factory is initialized with defaults.
%
% Last Saved Time-stamp: <Thu 2008-08-14 11:30:24 Eastern Daylight Time gramer>
%

  if ( ~isfield(stn, 'factories') || ~isfield(stn.factories, fac) )
    warning('No such factory yet! Starting with defaults...');
    stn.factories.(fac).variable = 'unspecified';
    stn.factories.(fac).fuzzies = ...
        { ...
            { 'drastic-low',	[ 0	1 ] }, ...
            { 'very-low',	[ 1	2 ] }, ...
            { 'low',		[ 2	3 ] }, ...
            { 'somewhat-low',	[ 3	4 ] }, ...
            { 'average',	[ 4	5 ] }, ...
            { 'somewhat-high',	[ 5	6 ] }, ...
            { 'high',		[ 6	7 ] }, ...
            { 'very-high',	[ 7	8 ] }, ...
            { 'drastic-high',	[ 8	9 ] }, ...
        };
  end;


  % Dialog options - see HELP INPUTDLG
  dlgopts.Interpreter = 'none';
  dlgopts.Resize = 'on';
  dlgopts.WindowStyle = 'modal';


  %
  % Get/change name of variable this factory will evaluate itself with
  %
  % (HACK: Prompt is long, so long var names will fit without resizing)
  variable = inputdlg('Variable (station field) that this factory will use to evaluate fuzzy values:', ...
                 'Variable name', ...
                 1, {stn.factories.(fac).variable}, dlgopts);
  % Cancel
  if ( isempty(variable) )
    warning('Edit cancelled...');
    return;
  end;
  % Validate variable name before proceeding
  if ( isempty(variable{:}) || ~iscellstr(variable) )
    error('What you entered was not a valid variable name string!');
  end;
  if ( ~isnan(str2double(variable{:}(1))) )
    error('Variable names must begin with a non-numeric character!');
  end;


  %
  % Get/change list of "fuzzy-value" strings factory might assign
  %
  % NOTE: Fuzzies 'unbelievably-low' and 'unbelievably-high' always ASSUMED
  fzs = [stn.factories.(fac).fuzzies{:}];
  fuzzarr = inputdlg('List of valid fuzzies (one per line, lowest to highest):', ...
                     'Fuzzy-Value Strings', ...
                     length(stn.factories.(fac).fuzzies), ...
                     {char(fzs{1:2:end})}, dlgopts);
  % Cancel
  if ( isempty(fuzzarr) )
    warning('Edit cancelled...');
    return;
  end;

  % Validate fuzzy strings before proceeding
  if ( isempty(strtrim(fuzzarr{:})) )
    error('Must have at least one fuzzy-value string!');
  end;

  fuzzies = cellstr(fuzzarr{:});
  if ( length(fuzzies) ~= length(unique(fuzzies)) )
    error('All fuzzy-value strings must be unique!');
  end;

  % Do we have the same number of fuzzies as we had before?
  endbnds = length(fzs);
  if ( size(fuzzies) ~= size(stn.factories.(fac).fuzzies) )
    endbnds = length(fuzzies) * 2;
    lastval = stn.factories.(fac).fuzzies{end}{2}(2);
    for idx = length(stn.factories.(fac).fuzzies)+1:length(fuzzies)
      fzs{(idx*2)-1} = fuzzies{idx};
      fzs{idx*2} = [lastval (lastval+1)];
      lastval = lastval + 1;
    end;
  end;


  %
  % Get/change numeric lower- and upper-bound values for each "fuzzy-value"
  %
  bnds = fzs(2:2:endbnds);
  bnds = cellstr(num2str(reshape([bnds{:}]', [2, length(bnds)])'));
  rngs = inputdlg(fzs(1:2:endbnds), ...
                  'Ranges', ...
                  1, bnds, dlgopts);
  % Cancel
  if ( isempty(rngs) )
    warning('Edit cancelled...');
    return;
  end;


  % Validate all ranges before proceeding
  for idx = 1:length(rngs)
    bds{idx} = sscanf(rngs{idx}, '%f %f');
    if ( length(bds{idx}) ~= 2 )
      error('Fuzzy %d (%s) needs two bounds - an upper and a lower!', ...
            idx, fuzzies{idx});
    end;
    if ( bds{idx}(1) >= bds{idx}(2) )
      error('Fuzzy %d (%s) bounds [%g,%g] are invalid!', ...
            idx, fuzzies{idx}, bds{idx}(1), bds{idx}(2));
    end;
    if ( idx > 1 && bds{idx}(1) ~= bds{idx-1}(2) )
      error('Fuzzy %d (%s) LOWER bound (%g) MUST equal previous (%s) UPPER (%g)!', ...
            idx, fuzzies{idx}, bds{idx}(1), fuzzies{idx-1}, bds{idx-1}(2));
    end;
  end;
  if ( isinf(bds{1}(1)) )
    warning('Lowest bound "-Inf" means NO value will ever be "unbelievably-low"!')
  end;
  if ( isinf(bds{end}(2)) )
    warning('Upper bound "+Inf" means NO value will ever be "unbelievably-high"!')
  end;
  sensrng = valid_var_range(variable{:});
  if ( sensrng(1) > bds{1}(1) || bds{end}(2) > sensrng(2) )
    warning('%s\n%s\n%s', ...
            sprintf('The range of valid values you specified was [%g,%g].', ...
                    bds{1}(1), bds{end}(2)), ...
            sprintf('\t"%s"', variable{:}), ...
            sprintf('ICON/matlab believes this variable is only valid between [%g,%g]!', ...
                    sensrng(1), sensrng(2)) ...
            );
  end;


  % Update stn.factories field for this factory with new data
  stn.factories.(fac).variable = variable{:};
  for idx = 1:length(rngs)
    stn.factories.(fac).fuzzies{idx} = { fuzzies{idx}, bds{idx}' };
  end;
  disp('Changes have been applied.');


  % Attempt to make sure we have some data for this (new? modified?) factory
  stn = verify_variable(stn, stn.factories.(fac).variable);
  if ( ~isfield(stn, stn.factories.(fac).variable) )
    warning('Could NOT create variable "%s" from existing data!', ...
            stn.factories.(fac).variable);
  end;

return;
