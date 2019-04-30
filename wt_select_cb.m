function output_txt = wt_select_cb(obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

pos = get(event_obj,'Position');
dvec = datevec(pos(1));
jday = num2str( fix(pos(1) - datenum(dvec(1), 1, 1) + 1) );
dat = 2.^(pos(2));
output_txt = {['X: ',datestr(pos(1)),' (',jday,')'],...
              ['Y: ',sprintf('%f', dat)]};

% %dat = 24*(10.^cast(pos(2),'double'))/(2*pi);
% %dat = 24*(10.^pos(2))/(2*pi);
% dat = (cast(pos(2),'double'))/(2*pi);
% output_txt = {['X: ',datestr(pos(1)),' (',jday,')'],...
%               ['Y: ',sprintf('%f', dat)]};

% % If there is a Z-coordinate in the position, display it as well
% if length(pos) > 2
%   output_txt{end+1} = ['Z: ',sprintf('%f', pos(3))];
% end

% % If there is CData, then also display it/them
% hdl = get(event_obj, 'Target');
% if ( ~isempty(hdl) )
%   xdata = get(hdl, 'XData');
%   ydata = get(hdl, 'YData');
%   zdata = get(hdl, 'ZData');
%   cdata = get(hdl, 'CData');
%   % If plot caller did not mesh X and Y, do so now
%   if (any(size(xdata) == 1))
%     [xdata,ydata] = meshgrid(xdata,ydata);
%   end;
%   if ( length(pos) >= 3 )
%     idx = find(xdata == pos(1) & ydata == pos(2) & zdata == pos(3),1);
%   else
%     idx = find(xdata == pos(1) & ydata == pos(2),1);
%   end; 
%   if ( isempty(idx) )
%     output_txt{end+1} = ['C: NO MATCH?'];
%   elseif ( isnumeric(cdata) && numel(cdata) >= idx )
%     output_txt{end+1} = ['C: ',sprintf('%f', cdata(idx))];
%   elseif ( isnumeric(cdata) && isscalar(cdata) )
%     output_txt{end+1} = ['C: ',sprintf('%f', cdata)];
%   elseif ( numel(zdata) >= idx )
%     output_txt{end+1} = ['C: ',sprintf('%f', zdata(idx))];
%   else
%     output_txt{end+1} = ['C: NO DATA?'];
%   end
% end
