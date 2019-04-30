function output_txt = xyz_map_select_cb(obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).
  pos = get(event_obj,'Position');
  [lat_loc,lon_loc] = inputm;
  output_txt = {['Lon: ',lon_loc],...
                ['Lat: ',lat_loc]};

  % If there is a Z-coordinate in the position, display it as well
  if length(pos) > 2
    output_txt{end+1} = ['Value: ',num2str(pos(3),4)];
  end;

return;
