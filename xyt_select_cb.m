function output_txt = xyt_select_cb(obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).
%
% Identical to TYZ_SELECT_CB, except this assumes that X and Y axis are
% simple numeric coordinates, and that Z is a vector of DATENUMs: this is
% useful, for example, in exploring interannual variability in data...

pos = get(event_obj,'Position');
dvec = datevec(pos(3));
jday = num2str( fix(pos(3) - datenum(dvec(1), 1, 1) + 1) );

pointno = '';
% If we can find A UNIQUE POINT corresponding to our selected data, then also
% show user the index in our XData vector corresponding to that point
eventattrs = get(event_obj);
if ( isfield(eventattrs,'Target') )
  tgt = eventattrs.Target;
  attrs=get(tgt);
  if ( all(isfield(attrs,{'XData','YData'})) ...
       && ~isempty(attrs.XData) && ~isempty(attrs.YData) ...
       && all(size(attrs.XData) == size(attrs.YData)) )
    pointix = find(pos(1) == attrs.XData & pos(2) == attrs.YData);
    if ( numel(pointix) == 1 )
      pointno = [' (X#' num2str(pointix) ')'];
    end;
  end;
end;

output_txt = {['X: ',sprintf('%g', pos(1))],...
              ['Y: ',sprintf('%g', pos(2)),pointno],...
              ['Z: ',datestr(pos(3)),' (',jday,')']...
             };
