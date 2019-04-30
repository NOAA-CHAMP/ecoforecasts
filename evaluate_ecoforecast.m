function events = evaluate_ecoforecast(ef, stn, begdate, enddate)
%function events = evaluate_ecoforecast(ef, stn, begdate, enddate)
%
% Evaluate ecological forecast 'ef' on station 'stn'. Returns a vector of
% 'event' structs, containing the set of datenums, conditions and fuzzy
% values matched by any rule(s) in the ecoforecast. If args 'begdate' and
% 'enddate' form a VALID daterange, only facts in that range are evaluated.
%
% EF is a struct with one field for each rule. Each rule is in turn a struct
% with two fields 'facts' and 'fuzzies': rule.facts a cell array of strings,
% each string naming a "fact-sensor"; rule.fuzzies is a cell array of cells,
% each element of rule.fuzzies being a list (cell array) of all acceptable
% "fuzzy" values for the corresponding fact-sensor in rule.facts. NOTE: If
% rule.fuzzies is an empty cell array, ANY fact of that type will match it,
% except the hard-coded values 'unbelievably-low' and 'unbelievably-high'.
%
% NOTE: For each string 'fact1' in rule.facts, there must exist as a field in
% the struct 'stn', a "fact-factory" struct named 'stnf1.factories.fact1',
% which can be used to convert data from some time series that is a field in
% that station, into a cell array of symbolic/semantic values ("facts").
% This is known as "factizing" a variable (q. EVALUATE_FACTORY).
%
% Last Saved Time-stamp: <Tue 2010-02-16 19:03:45 Eastern Standard Time gramer>


  events.name = ef.name;
  events.shortdesc = ef.shortdesc;
  events.desc = ef.desc;

  events.jday = [];
  alljdays = [];
  allsris = [];


  if ( ~exist('begdate', 'var') || ~exist('enddate', 'var') || ...
       isempty(begdate) || isempty(enddate) || begdate >= enddate )
    begdate = [];
    enddate = [];
  end;

  if ( ~isstruct(ef) || ~isstruct(stn) )
    error('Both arguments should be structs!');
  end;
  if ( ~isfield(ef, 'rules') )
    error('No rules have been specified in "ef"!');
  end;
  if ( ~isfield(stn, 'factories') )
    error('No fact factories have been specified in "stn"!');
  end;
  if ( length(fieldnames(stn)) < 2 )
    error('Struct "stn" has no sensor data (only one field)!');
  end;


  % Sort rules of EF based on "salience" (most complex==most facts, first)
  rules = fieldnames(ef.rules);
  for ri = 1:length(rules)
    rule = ef.rules.(rules{ri});
    % Do not count "catch-all" rules - those with empty fuzzy lists!
    salience = sum(~cellfun(@isempty,ef.rules.(rules{ri}).fuzzies));
    rulearr{ri} = sprintf('%03d %s', salience, rules{ri});
  end;
  rules = cellstr(sortrows(char(rulearr), [-1 -2 -3]));
  for ri = 1:length(rules)
    rules{ri} = rules{ri}(5:end);
  end;


  % Construct our 'agenda' of all facts that are relevant to any rule
  % Note 'agenda' will have only one field per fact type ('fact sensor'),
  % no matter how many of our rules may rely on that particular fact type.
  agenda = [];
  for ri = 1:length(rules)
    rule = ef.rules.(rules{ri});
    if ( ~isfield(rule, 'facts') )
      error('No fact sensors in rule (no ef.rules.%s.facts)!', rules{ri});
    end;
    if ( ~isfield(rule, 'fuzzies') )
      error('No fuzzies in rule (no ef.rules.%s.fuzzies)!', rules{ri});
    end;

    % Allow for conditions that may match one of several 'similar'
    % facts, e.g., a sea temperature condition that can match fact
    % sensors named either 'seandbc' and/or 'sea1m', etc. I.e., if
    % an element of cell-array rule.facts is itself a cell-array of
    % strings, then handle each string in that array in turn...
    for fi = 1:length(rule.facts)

      facts = rule.facts{fi};
      if ( ~iscell(facts) )
        facts = { facts };
      end;
      foundAnyFactory = false;
      for fii = 1:length(facts)
        fact = facts{fii};
        if ( isfield(stn.factories, fact) )
          foundAnyFactory = true;
          if ( isempty(agenda) || ~isfield(agenda, fact) )
            agenda.(fact) = evaluate_factory(stn, fact, [], ...
                                             begdate, enddate);
          end;
        else
          agenda.(fact).datenum = [];
          agenda.(fact).averages = [];
          agenda.(fact).fuzzies = {};
        end;
      end;
      if ( ~foundAnyFactory )
        warning('No factory in stn.factories for rule "%s", fact(s): "%s"!', ...
                rules{ri}, sprintf('%s,', facts{:}));
      end;
    end;
  end;

  %%%% DEBUG:  save('DEBUG_agenda.mat', 'agenda');

  % Match facts from our agenda with rules in salience order, then remove
  % each matching fact so that lower-precedence rules will not see it
  for ri = 1:length(rules)

    rule = ef.rules.(rules{ri});

    clear hits;
    hits(1).datenum = [];
    hits(1).fuzzies = {};
    hits(1).durations = [];
    hits(1).sri = [];
    hits(2:length(rule.facts)) = hits(1);

    % For this rule, evaluate each of its conditions in turn
    % NOTE: Salience is not meaningful here: if two conditions are both
    % competing for the same kind of fact, evaluation order is random.
    curhit = 0;
    for fi = 1:length(rule.facts)

      facts = rule.facts{fi};
      fuzzies = rule.fuzzies{fi};

      if ( ~iscell(facts) )
        facts = { facts };
      end;
      for fii = 1:length(facts)

        fact = facts{fii};
        curhit = curhit + 1;

        % CATCH-ALL: If rule specifies no fuzzies for this type of fact, then
        % accept ALL of them, except the unbelievable values.
        if ( isempty(fuzzies) )
          hitidx = find(~ismember(agenda.(fact).fuzzies, ...
                                  {'unbelievably-low','unbelievably-high'}));
        else
          hitidx = find(ismember(agenda.(fact).fuzzies, fuzzies));
        end;
        hits(curhit).datenum = agenda.(fact).datenum(hitidx);
        hits(curhit).fuzzies = agenda.(fact).fuzzies(hitidx);
        if ( fii == 1 )
          factjdays = unique(floor(hits(curhit).datenum));
        else
          factjdays = union(factjdays, floor(hits(curhit).datenum));
        end;

      end;

      if ( fi == 1 )
        rulehits(ri).jday = factjdays;
      else
        rulehits(ri).jday = intersect(rulehits(ri).jday, factjdays);
      end;

      % If ANY combination of conditions is not met, rule FAILS
      if ( isempty(rulehits(ri).jday) )
        break;
      end;

    end;


    % Were there "events" - i.e., any days when ALL conditions matched?
    if ( ~isempty(rulehits(ri).jday) )

      curhit = 0;
      for fi = 1:length(rule.facts)

        facts = rule.facts{fi};

        if ( ~iscell(facts) )
          facts = { facts };
        end;
        for fii = 1:length(facts)

          fact = facts{fii};
          curhit = curhit + 1;

          % Keep only hits from jdays when ALL conditions matched
          hitidx = find(ismember(floor(hits(curhit).datenum), rulehits(ri).jday));

          hits(curhit).datenum = hits(curhit).datenum(hitidx);
          hits(curhit).fuzzies = hits(curhit).fuzzies(hitidx);

          if ( isempty(hitidx) )

            % If this fact-sensor (of this condition) had no matches...
            events.rules.(rules{ri}).(fact).datenum = [];
            events.rules.(rules{ri}).(fact).fuzzies = [];
            events.rules.(rules{ri}).(fact).durations = [];
            events.rules.(rules{ri}).(fact).sri = [];

          else

            % Take only the longest contiguous GROUP of matches each day;
            % leave the rest available for other rules to match against.
            kpidx = [];
            rmidx = [];

            % First break up the vector of hits into individual days
            % ??? UT/GMT: This will be WRONG for Indo-Pacific sites!
            dayix = find(diff(floor(hits(curhit).datenum)) ~= 0);
            edays = [ 0 dayix(:)' length(hits(curhit).datenum) ] + 1;
            for di = 1:length(edays)-1
              %dn = hits(curhit).datenum(edays(di):edays(di+1)-1);
              fz = hits(curhit).fuzzies(edays(di):edays(di+1)-1);

              % Get length of contiguous indices having the same fuzzy
              % (Run-length encoding algorithm courtesy "MATLAB array
              % manipulation tips and tricks", Peter Acklam, Norway)
              kern = find(~strcmp(fz(1:end-1), fz(2:end)));
              rl = diff([ 0 kern(:)' length(fz) ]);
              ez = fz([ kern(:)' length(fz) ]);

              % Find the longest run-length: that's our "super-period"!
              [maxrl, maxrli] = max(rl);
              % Which indices in original hits datenum vector was that?
              begi = edays(di) + sum(rl(1:maxrli-1));
              contigi = begi:(begi + maxrl - 1);

              % Assign the S/RI to the FIRST matching fact in the group
              % (S/RI = Stimulus/Response Index for this ecoforecast)
              hits(curhit).durations(contigi(1),1) = maxrl * (3/24);
              hits(curhit).sri(contigi(1),1) = maxrl * fuzzy_sri(ez{maxrli});

              % Accumulate all the periods/super-periods we've chosen
              kpidx = [kpidx contigi(1)];
              rmidx = [rmidx contigi];
            end;

            % Remove every "super-fact" (i.e., including any other individual
            % facts accumulated into it!) that "fully matches" (i.e., when
            % all conditions from this rule were matched from the agenda
            agrmidx = find(ismember(agenda.(fact).datenum, hits(curhit).datenum(rmidx)));
            agenda.(fact).datenum(agrmidx) = [];
            agenda.(fact).fuzzies(agrmidx) = [];
            agenda.(fact).averages(agrmidx) = [];

            % Release all the other facts for following rules to use
            hits(curhit).datenum = hits(curhit).datenum(kpidx);
            hits(curhit).fuzzies = hits(curhit).fuzzies(kpidx);
            hits(curhit).durations = hits(curhit).durations(kpidx);
            hits(curhit).sri = hits(curhit).sri(kpidx);

            % Each matching fact we've kept, store inside the event...
            events.rules.(rules{ri}).(fact) = hits(curhit);

          end;

        end;

      end;

      % Accumulate all days (floored datenums) and S/RIs for this rule
      [xc,yc] = consolidator(floor([hits(:).datenum]), [hits(:).sri], 'sum');
      events.rules.(rules{ri}).jday = xc';
      events.rules.(rules{ri}).sri = yc';

      fprintf(1, '\t EF.%s: %d matching days\n', ...
              upper(rules{ri}), length(events.rules.(rules{ri}).jday));

      alljdays = [ alljdays events.rules.(rules{ri}).jday ];
      allsris = [ allsris events.rules.(rules{ri}).sri ];

    end;

  end;

  if ( ~isempty(alljdays) )
    [xc, yc] = consolidator(alljdays, allsris, 'sum');
    events.jday = xc';
    events.sri = yc';
  end;

  % Finally, calculate 'periods' of multiple subsequent event-days
  % (A 'period' is a sequence of two or more days, contiguous or nearly so
  % with each other, on which events of the same class occurred. And what
  % is or is not 'nearly contiguous' may depend on the ecoforecast...)
  if ( ~isempty(events.jday) )
    jdix = find( diff(events.jday) > 2 );
    idx = [ 1 ( jdix(:)' + 1 ) ];
    events.period = events.jday(idx);
  end;


return;
