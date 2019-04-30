function wpts = read_waypoint_xml(fname)
%function wpts = read_waypoint_xml(fname)
%
% Read GPS waypoints from the ".gpx" (XML-format) file FNAME. If available,
% also reads optional REEF_NAME from each waypoint node, as: "<extensions>"
% -> "<gpxx:WaypointExtension>" -> "<gpxx:Categories>" -> "<gpxx:Category>".
%
% Returns STRUCT vector WPTS, containing fields .station_name, .lon, .lat,
% (DATENUM) .timestamp, (numeric) .waypoint_num, and (optional) .reef_name.
%
% CALLS: XML2STRUCT (BioInfo Toolbox)
%
% Last Saved Time-stamp: <Fri 2018-10-26 13:44:50 Eastern Daylight Time gramer>

  str = xml2struct(fname);
  
  cnms = string({str.Children.Name});
  wptix = find(cnms == 'wpt');

  wptstr = str.Children(wptix);
  for ix=1:numel(wptstr);
    wpts(ix).station_name = '';
    attrstrs = string({wptstr(ix).Attributes.Name});
    lonix = find(attrstrs == 'lon');
    latix = find(attrstrs == 'lat');
    if ( isempty(lonix) || isempty(latix) )
      error('Missing LON or LAT Attribute: WPT #%d, file %s',ix,fname);
    end;
    wpts(ix).lon = str2num(wptstr(ix).Attributes(lonix).Value);
    wpts(ix).lat = str2num(wptstr(ix).Attributes(latix).Value);

    chldstr = string({wptstr(ix).Children.Name});
    namix = find(chldstr == 'name');
    timix = find(chldstr == 'time');
    if ( isempty(namix) || isempty(timix) )
      error('Missing NAME or TIME Attribute: WPT #%d, file %s',ix,fname);
    end;
    wpts(ix).station_name = string(wptstr(ix).Children(namix).Children(1).Data);
    wpts(ix).waypoint_num = str2num(wptstr(ix).Children(namix).Children(1).Data);
    wpts(ix).timestamp = datenum(wptstr(ix).Children(timix).Children(1).Data,'yyyy-mm-ddTHH:MM:SSZ');

    extix = find(chldstr == 'extensions');
    if ( ~isempty(extix) )
      extnms = string({wptstr(ix).Children(extix).Children.Name});
      gwxix = find(extnms == 'gpxx:WaypointExtension');
      if ( ~isempty(gwxix) )
        extstr = string({wptstr(ix).Children(extix).Children(gwxix).Children.Name});
        rfnix = find(extstr == 'gpxx:Categories');
        if ( ~isempty(rfnix) )
          wpts(ix).reef_name = ...
              string(wptstr(ix).Children(extix).Children(gwxix).Children(rfnix).Children(2).Children(1).Data);
        end;
      end;
    end;
  end;

return;
