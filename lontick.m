function lontick(ah, axnm, selfun)
% LONTICK  Tick labels with negative longitudes in place of (180,360]
%
%function lontick(ah, axnm, selfun)

  if (nargin == 1); if (~ishandle(ah)); axnm=ah; ah=[]; end; end;

  if ( ~exist('ah','var') || isempty(ah) ); ah = gca; end;
  if ( ~exist('axnm','var') || isempty(axnm) ); axnm = 'x'; end;

  if ( ~exist('selfun','var') || isempty(selfun) );
    switch (lower(axnm))
     case 'x',      selfun = @lonyzc_select_cb;
     case 'y',      selfun = @xlonzc_select_cb;
     otherwise,     error('Axis name must be "x" or "y"!');
    end;
  end;

  % Required for pan/zoom to work properly??? (Q.v. DATETICK 'keeplimits', 'keepticks')
  set(ah, [axnm 'tickmode'],'auto');
  set(ah, [axnm 'ticklabelmode'], 'auto');

  % Initialize callback data struct (borrowed liberally from DATETICK2, q.v.)
  ltd.ah = ah;
  ltd.axnm = axnm;
  ltd.selfun = selfun;

  % Get the handles to other parts of the figure
  ltd.parent=get(ltd.ah,'parent');  % The figure handle
  ltd.kids=get(ltd.parent,'Children');  % All children of the figure
  ltd.axes=[];   % Children of the figure that are also axes
  for n=1:length(ltd.kids)
    if strcmp(get(ltd.kids(n),'type'),'axes') ...
          && ~strcmp(get(ltd.kids(n),'tag'),'Legend') ...
          && ~strcmp(get(ltd.kids(n),'tag'),'Colorbar')
      ltd.axes=[ltd.axes ltd.kids(n)];
    end
  end
  % Link all the axes together
  ltd.hlink=linkprop(ltd.axes,{[ltd.axnm 'lim'];[ltd.axnm 'tick']});
  ltd.hlink2=linkprop(ltd.axes,[ltd.axnm 'ticklabel']); %for some reason this can't be on the line above

  % Set application data and callbacks so ticks will update after pan and zoom
  for n=1:length(ltd.axes)
    setappdata(ltd.axes(n),'lontickdata',ltd);
  end
  zh = zoom(ltd.parent);
  if (ishandle(zh)); set(zh,'ActionPostCallback',@lontick_hdlr); end;
  ph = pan(ltd.parent);
  if (ishandle(ph)); set(ph,'ActionPostCallback', @lontick_hdlr); end;
  % If user clicks a point while in Data Cursor mode, report "right" longitude
  dh = datacursormode(ltd.parent);
  if (ishandle(dh)); set(dh,'UpdateFcn',selfun); end;


  % Do actual tick label modification
  lons = get(ah, [axnm 'tick']);

  meridix = find(lons > 180);
  lons(meridix) = lons(meridix) - 360;
  lonlabels = num2str(lons(:));

  % set(ah, [axnm 'tick'], lons, [axnm 'ticklabel'], lonlabels);
  set(ah, [axnm 'ticklabel'], lonlabels);

return;



%%%%%%%%%% 
%%%%%%%%%% PRIVATE FUNCTIONS (CALLBACKS)
%%%%%%%%%% 


function lontick_hdlr(varargin)
%function lontick_hdlr(varargin)
% Reformat longitude axis (>180 => negative) after Pan or Zoom

    ltd = getappdata(varargin{2}.Axes,'lontickdata');
    ah = gca; %used to be ltd.axh instead of gca. But after a zoom, we always want the current axes

    % Redo lontick(ah, ltd.axnm, ltd.selfun) without setting callbacks again...
    % I.e., do actual tick label modification
    lons = get(ah, [ltd.axnm 'tick']);

    meridix = find(lons > 180);
    lons(meridix) = lons(meridix) - 360;
    lonlabels = num2str(lons(:));

    % set(ah, [axnm 'tick'], lons, [axnm 'ticklabel'], lonlabels);
    set(ah, [ltd.axnm 'ticklabel'], lonlabels);
    set(ah, [ltd.axnm 'tickmode'], 'auto');
    set(ah, [ltd.axnm 'ticklabelmode'], 'auto');

return;




function output_txt = lonyzc_select_cb(obj,event_obj)
%function output_txt = lonyzc_select_cb(obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

  pos = get(event_obj,'Position');

  lon = pos(1);
  if (lon > 180); lon = lon - 360; end;

  output_txt = {['X: ',num2str(lon)],...
                ['Y: ',num2str(pos(2))]};

  output_txt = append_ZC(output_txt,obj,event_obj,pos);

return;




function output_txt = xlonzc_select_cb(obj,event_obj)
%function output_txt = xlonzc_select_cb(obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

  pos = get(event_obj,'Position');

  lon = pos(2);
  if (lon > 180); lon = lon - 360; end;

  output_txt = {['X: ',num2str(pos(1))],...
                ['Y: ',num2str(lon)]};

  output_txt = append_ZC(output_txt,obj,event_obj,pos);

return;


function output_txt = append_ZC(output_txt,obj,event_obj,pos)
%function output_txt = append_ZC(output_txt,obj,event_obj,pos)

  % If there is a Z-coordinate in the position, display it as well
  if numel(pos) > 2
    output_txt{3} = ['Z: ',sprintf('%f', pos(3))];
  end

  % Try to also display X/Y (and Z) INDICES as well as values...
  hdl = get(event_obj, 'Target');
  % if ( ~isempty(hdl) )
  if ( ishandle(hdl) )
      xdata = get(hdl, 'XData');
      ydata = get(hdl, 'YData');
      zdata = get(hdl, 'ZData');
      try,
        cdata = get(hdl, 'CData');
      catch,
        cdata = '';
      end;

      idx = [];

      if ( isempty(zdata) || numel(pos) < 3 )
        % If plot caller did not mesh X and Y, do so now
        if (any(size(xdata) == 1))
          [xdata,ydata] = meshgrid(xdata,ydata);
        end;
        idx = find(xdata == pos(1) & ydata == pos(2));
        [ix,jx] = ind2sub(size(xdata), idx);
        output_txt{1} = sprintf('%s \t (col %d)', output_txt{1}, ix);
        output_txt{2} = sprintf('%s \t (row %d)', output_txt{2}, jx);
      else
        % If plot caller did not mesh X/Y (with or w/o Z), do so now
        if (any(size(xdata) == 1))
          if (any(size(zdata) == 1))
            [xdata,ydata,zdata] = meshgrid(xdata,ydata,zdata);
          else
            [xdata,ydata] = meshgrid(xdata,ydata);
          end;
        end;
        idx = find(xdata == pos(1) & ydata == pos(2) & zdata == pos(3));
        [ix,jx,kx] = ind2sub(size(xdata), idx);
        output_txt{1} = sprintf('%s \t (col %d)', output_txt{1}, ix);
        output_txt{2} = sprintf('%s \t (row %d)', output_txt{2}, jx);
        if ( kx ~= 1 )
          output_txt{3} = sprintf('%s \t (lyr %d)', output_txt{3}, kx);
        end;
      end;

      % If there is CData (e.g., this graph was produced by scatter or scatter3,
      % instead of plot or plot3), then find the C data and also display it/them.
      if ( isempty(idx) || numel(cdata) < idx )
        % output_txt{4} = ['C: NO MATCH?'];
      else
        % Note there may be MORE THAN ONE matching C value output by SPRINTF
        output_txt{4} = ['C: ',sprintf('%f ', cdata(idx))];
      end
  end

return;
