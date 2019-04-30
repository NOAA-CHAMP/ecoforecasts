1;

doPrint = true;
%doPrint = false;

methods = {'mellor','monismith','shay','ardhuin'};
%waveper = [3,4,6,8,10];
waveper = [3,6,9];

ws = linspace(  2, 20,20);
hs = linspace(0.3,1.8,20);
[WS,HS] = meshgrid(ws,hs); 

% Stokes drift in kilometers per day
STkpd = repmat(nan,[numel(methods),numel(waveper),size(WS,1),size(WS,2)]);

for mthix=1:numel(methods);
  method = methods{mthix};

  fmg;
  for wvix=1:numel(waveper);
    per = waveper(wvix);
    %ST = stokes_drift(WS,HS,per,method);
    ST = stokes_drift_2(WS,HS,per,method);
    z = ST*3600*24/1e3;
    surf(WS,HS,z);
    text(WS(end,end),HS(end,end),z(end,end),['P=',num2str(per),' s']);
    % Save it for later / Don't run away and let me down --English Beat
    STkpd(mthix,wvix,:,:) = z;
  end;
  view(45,45);
  xlabel('Wind speed [kts]');
  ylabel('Sig. wave hgt. [m]');
  zlabel('Current speed [km^.d^-^1]');
  titlename(['Surface current - method of ',capitalize(method),' et al. (2)']);
  if doPrint; print('-dpng',[mfilename,'_2_',method,'.png']); end;
end;
