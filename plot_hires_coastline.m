function [ch,ax] = plot_hires_coastline(bath,ax,useHighestRes)
%function [ch,ax] = plot_hires_coastline(bath,ax,useHighestRes)
%
% Plot double-wide black line for coastline within region covered by the
% bathymetry STRUCT BATH: BATH should have fields .lon,.lat,.field [m]; it
% may also be a cell array of matrices {LON,LAT,FIELD}.
%
% Plots in the current axes (GCA), unless optional arg AX is specified.
%
% If either USEHIGHESTRES (DEFAULT: False) or field is entirely within
% south Florida, load file SOFLA_COAST.dat and call FILL to contour it. 
% Otherwise, CONTOUR [0,0], or if no 0-contour is found, CONTOUR [-1,-1].
%
% Optionally returns handle CH returned by FILL or CONTOUR, and AXES AX.
%
% Last Saved Time-stamp: <Sat 2017-04-22 12:56:52 Eastern Daylight Time gramer>

  global sofla_coast;

  if ( ~exist('bath','var') )
    bath=[]; % THIS IS AN ERROR - we should catch it below
  end;
  if ( iscell(bath) && numel(bath)==3 && ...
       ismatrix(bath{1}) && ismatrix(bath{2}) && ismatrix(bath{3}) )
    cel = bath;
    bath=[]; clear bath
    bath.lon = cel{1};
    bath.lat = cel{2};
    bath.field = cel{3};
    cel=[]; clear cel
  end;
  if ( isempty(bath) )
    error('First arg BATH must be a bathymetry STRUCT or 3-element cell array');
  end;

  if ( ~exist('ax','var') || isempty(ax) )
    ax = gca;
  end;
  if ( ~exist('useHighestRes','var') || isempty(useHighestRes) )
    useHighestRes = false;
  end;

  bbox = [min(bath.lon(:)),max(bath.lon(:)),...
          min(bath.lat(:)),max(bath.lat(:))];
  % For locations in south Florida, use high-resolution coastline
  % Elsewhere, make do with the high-resolution 0 m isobath
  if ( ~useHighestRes || bboxint(bbox,[-82.9261,-80.0310,24.5017,27.0000]) < eps )
    [cch,ch] = contour(ax,bath.lon,bath.lat,...
                      bath.field,[0,0],'k-','LineWidth',2);
    if ( isempty(cch) )
      % Some high-resolution bathymetry datasets don't QUITE reach the coast:
      % as a last-ditch effort, use the -1 m isobath in place of 0.
      [cch,ch] = contour(ax,bath.lon,bath.lat,...
                         bath.field,[-1,-1],'k-','LineWidth',2);
    end;
  else
    if ( ~exist('sofla_coast', 'var') || isempty(sofla_coast) )
      disp('Reloading coastline');
      % sofla_coast = load('sofla_coast_low.dat');
      % sofla_coast = load('sofla_coast_medium.dat');
      sofla_coast = load('sofla_coast.dat');
    end;
    axes(ax);
    %%fill(sofla_coast(:,1), sofla_coast(:,2), [.4 .3 .2]);
    %ch=fill(sofla_coast(:,1), sofla_coast(:,2), [0.0 0.0 0.0], 'LineWidth',2);
    useix = find(bbox(1)<=sofla_coast(:,1)&sofla_coast(:,1)<=bbox(2) & bbox(3)<=sofla_coast(:,2)&sofla_coast(:,2)<=bbox(4));
    ch=fill(sofla_coast(useix,1), sofla_coast(useix,2), [0.0 0.0 0.0], 'LineWidth',2);
    axis(bbox);
  end;

  daspect([1,cosd(bath.lat(1)),1]);

  % In datamode, user is interested in bathymetry or other contours, not coastline
  set(ch,'HitTest','off')

return;
