function sri = fuzzy_sri(fuzzy)

% PER ICON/G2
  sri = 3;
  switch ( fuzzy ),
   case {'drastic-low', 'drastic-high', 'conducive'},
    sri = sri * 2.5;
   case {'very-low', 'very-high', 'somewhat-conducive'},
    sri = sri * 2;
  end;

% % PER Jim's (and Jank's) original CLIPS code
%   sri = 3;
%   switch ( fuzzy ),
%    case {'drastic-low', 'drastic-high', 'conducive'},
%     sri = sri * 2;
%   end;

return;
