1;

if ( ~exist('doPrint','var') || isempty(doPrint) )
  doPrint = false;
end;
if ( ~exist('upwpath','var') )
  upwpath = fullfile(get_coral_path,'CRCP','Upwelling','CoRIS');
end;


for dj=0:7:56;
  disp(datestr(xl(1)+dj));
  figbasenm = fullfile(upwpath,[mfilename,'_',datestr(xl(1)+dj,'yyyymmmdd')]);

  if ishandle(tpax); axes(tpax(end)); xlim(xl+dj); datetick3;
    if doPrint; print('-dpng',[figbasenm,'_temps.png']); end; end;
  pause(0.5);
  if ishandle(xsax); axes(xsax(end)); xlim(xl+dj); datetick3;
    if doPrint; print('-dpng',[figbasenm,'_cross.png']); end; end;
  pause(0.5);
  if ishandle(lsax); axes(lsax(end)); xlim(xl+dj); datetick3;
    if doPrint; print('-dpng',[figbasenm,'_along.png']); end; end;
  pause(0.5);
  if ishandle(dpax); axes(dpax(end)); xlim(xl+dj); datetick3;
    if doPrint; print('-dpng',[figbasenm,'_deep.png']); end;
  pause(0.5);
  if ishandle(shax); axes(shax(end)); xlim(xl+dj); datetick3;
    if doPrint; print('-dpng',[figbasenm,'_shallow.png']); end; end;
  pause(0.5);
  if ishandle(ekax); axes(ekax(end)); xlim(xl+dj); datetick3;
    if doPrint; print('-dpng',[figbasenm,'_ekman.png']); end; end;
  pause(0.5);
  if ishandle(wsax); axes(wsax(end)); xlim(xl+dj); datetick3;
    if doPrint; print('-dpng',[figbasenm,'_stress.png']); end; end;
  pause(0.5);
  keyboard;
end;
