function mdl = insolation_model(t, site, unit_conv)
%function mdl = insolation_model(t, site, unit_conv)
%
% Implement van Woesik et al (2006) two-sine curve insolation model
%
% Last Saved Time-stamp: <Wed 2008-08-13 14:45:24 Eastern Daylight Time gramer>
%

    a1 = site(1); b1 = site(2); c1 = site(3);
    a2 = site(4); b2 = site(5); c2 = site(6);

    mdl = ( a1 * sin((b1 * t) + c1) ) ...
        + ( a2 * sin((b2 * t) + c2) );

    if ( exist('unit_conv', 'var') && ~isempty(unit_conv) && isnumeric(unit_conv) )
      mdl = mdl * unit_conv;
    end;

return;
