function [variables, sensors] = associate_variables(stn, events)
%function [variables, sensors] = associate_variables(stn, events)
%
% Find the names of all station data variables associated with ecoforecast
% events in struct 'events' (normally, returned by evaluate_ecoforecast).
% Also returns cell array of same size as 'variables', containing the
% 'fact-sensor' names of all matched facts in the events list.
%
% Last Saved Time-stamp: <Thu 2009-01-22 09:45:05 Eastern Standard Time gramer>
%

    variables = {};

    rules = fieldnames(events.rules);

    sensors = {};
    for ri = 1:length(rules)
        sensors = union(sensors, fieldnames(events.rules.(rules{ri})));
    end;

    sensors = unique(sensors);
    sensors(strcmp(sensors, 'jday')) = [];

    vari = 0;
    vars = {};
    rmidx = [];
    facs = fieldnames(stn.factories);
    for si = 1:length(sensors)
        findx = find( strcmp(facs, sensors{si}), 1 );
        if ( isempty(findx) )
            rmidx = [ rmidx findx ];
        else
            vari = vari + 1;
            vars{vari} = stn.factories.(facs{findx}).variable;
            sens{vari} = sensors{si};
        end;
    end;

    [variables, vari] = unique(vars);
    sensors = sens(vari);

return;
