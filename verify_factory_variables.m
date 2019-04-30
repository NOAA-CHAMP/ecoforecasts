function stn = verify_factory_variables(stn)
%function stn = verify_factory_variables(stn)
%
% Make sure that each field in struct 'stn.factories' names a valid variable,
% i.e., one that has been loaded or calculated already. Calls VERIFY_VARIABLE.
%
%
% Last Saved Time-stamp: <Thu 2008-08-14 07:13:23 Eastern Daylight Time gramer>
%

    if ( ~isfield(stn, 'factories') )
        warning('Struct "stn" has no field "stn.factories"!');
    else
        sensors = fieldnames(stn.factories);
        for idx = 1:length(sensors)
            sensor = sensors{idx};
            variable = stn.factories.(sensor).variable;
            stn = verify_variable(stn, variable);
        end;
    end;

return;
