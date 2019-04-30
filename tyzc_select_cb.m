function output_txt = tyzc_select_cb(obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

pos = get(event_obj,'Position');
dvec = datevec(pos(1));
jday = num2str( fix(pos(1) - datenum(dvec(1), 1, 1) + 1) );
output_txt = {['X: ',datestr(pos(1)),' (',jday,')'],...
    ['Y: ',sprintf('%f', pos(2))]};

% If there is a Z-coordinate in the position, display it as well
if length(pos) > 2
    output_txt{end+1} = ['Z: ',sprintf('%f', pos(3))];
end

% If there is CDate (e.g., this graph was produced by scatter or scatter3,
% instead of plot or plot3), then find the C data and also display it/them.
hdl = get(event_obj, 'Target');
if ( ~isempty(hdl) )
    xdata = get(hdl, 'XData');
    ydata = get(hdl, 'YData');
    zdata = get(hdl, 'ZData');
    try,
      cdata = get(hdl, 'CData');
    catch,
      cdata = '';
    end;
    % If plot caller did not mesh X and Y, do so now
    if (any(size(xdata) == 1))
      [xdata,ydata] = meshgrid(xdata,ydata);
    end;
    if ( length(pos) >= 3 )
      idx = find(xdata == pos(1) & ydata == pos(2) & zdata == pos(3),1);
    else
      idx = find(xdata == pos(1) & ydata == pos(2),1);
    end; 
    if ( isempty(idx) )
        output_txt{end+1} = ['C: NO MATCH?'];
    elseif ( ~ischar(cdata) && numel(cdata) >= idx )
        output_txt{end+1} = ['C: ',sprintf('%f', cdata(idx))];
    elseif ( numel(zdata) >= idx )
        output_txt{end+1} = ['C: ',sprintf('%f', zdata(idx))];
    else
        output_txt{end+1} = ['C: NO DATA?'];
    end
end
