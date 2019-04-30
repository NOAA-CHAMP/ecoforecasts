function hits = eval_ef(ef, stn)
%function hits = eval_ef(ef, stn)
%
% Evaluate ecological forecast 'ef' on station 'stn'. Returns vector 'hits',
% the set of all datenums on which a match for any rule was found. Note that
% size(hits) == size(fieldnames(ef)) always.
%
% EF is a struct with one field for each rule. Each rule is in turn a struct
% with two fields 'facts' and 'fuzzies': rule.facts a cell array of strings,
% each string naming a "fact-sensor"; rule.fuzzies is a cell array of cells,
% each element of rule.fuzzies being a list (cell array) of all acceptable
% "fuzzy" values for the corresponding fact-sensor in rule.facts.
%
% NOTE: For each string 'fact1' in rule.facts, there must exist as a field in
% the struct 'stn', a "fact-factory" struct named 'stnf1.factories.fact1',
% which can be used to convert data from some time series that is a field in
% that station, into a cell array of symbolic values ("facts"). The result is
% a series of strings, one for each value in the time-series.
%
% Last modified: Lew.Gramer@noaa.gov, 2008 Jul 31
%

    if ( ~isfield(stn, 'factories') )
        error('No fact factories have been specified in "stn"!');
    end;


    rules = fieldnames(ef);

    for ri = 1:length(rules)

        rule = ef.(rules{ri});
        hits(ri) = [];

        if ( ~isfield(rule, 'facts') )
            error('No facts found in rule ef.%s!', rules{ri});
        end;
        if ( ~isfield(rule, 'fuzzies') )
            error('No fuzzies found in rule ef.%s!', rules{ri});
        end;


        for fi = 1:length(rule.facts)

            fact = rule.facts{fi};
            fuzzies = rule.fuzzies{fi};

            if ( ~isfield(stn.factories, fact) )
                error('No factory in stn.factories for fact "%s"!', fact);
            end

            facthits = find(stn.factories.(fact)

            var = stn.(ef.(rules{ri}).variable){vi};
            lbd = ef.(rules{ri}).lower{vi};
            ubd = ef.(rules{ri}).upper{vi};
            if ( isempty(hits(ri)) )
              hits(ri) = find(lbd <= var & var <= ubd));
            else
              hits(ri) = intersect(hits(ri), ...
                                   find(lbd <= var & var <= ubd));
            end;

        end;

    end;

return;

%function dnum = find_matching_dates(dts, ts, 
