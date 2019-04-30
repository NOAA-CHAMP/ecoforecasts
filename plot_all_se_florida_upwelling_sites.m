1;
% SCRIPT to plot LON and LAT of all station and mooring sites whose data were
% analyzed in the NOAA AOML CHAMP Upwelling project (CRCP#789). Assumes the
% current Figure and AXES contain a map suitable for PLOT(lon,lat...) calls.
% Calls ANNOTAXIS (v.) to label site locations to the *left* of the map.
%
% Last Saved Time-stamp: <Tue 2018-04-10 13:45:12 EDT lew.gramer>

    if ( ~exist('fontsz','var') || isempty(fontsz) )
      %fontsz = 16;
      fontsz = 24;
    end;
    if ( ~exist('blowupmap','var') || isempty(blowupmap) )
      blowupmap = false;
    end;


    plotargs = {'rp','MarkerFaceColor','w','MarkerSize',10};

    %annotargs = {'FontSize',fontsz,'BackgroundColor','w'};
    annotargs = {'FontSize',fontsz,'BackgroundColor','none'};

    % Fix the AXES bounds, so stations outside the mapped area are ignored
    axis(axis);

    if ( exist('fwyf1','var') )
      plot(fwyf1.lon,fwyf1.lat,plotargs{:});
      annotaxis(fwyf1.lat,'FWYF1: met. \rightarrow',annotargs{:});
    end;
    if ( exist('pvgf1','var') )
      plot(pvgf1.lon,pvgf1.lat,plotargs{:});
      if ( blowupmap )
        annotaxis(pvgf1.lat,'Port Everglades: met. \rightarrow',annotargs{:});
      else
        annotaxis(pvgf1.lat+0.03,'Port Everglades: met. \rightarrow',annotargs{:});
      end;
    end;
    if ( exist('lkwf1','var') )
      plot(lkwf1.lon,lkwf1.lat,plotargs{:});
      annotaxis(lkwf1.lat,'LKWF1: met. \rightarrow',annotargs{:});
    end;
    if ( exist('fdepk','var') )
      plot(fdepk.lon,fdepk.lat,plotargs{:});
      %annotaxis(fdepk.lat,'FDEP "K": met. \rightarrow',annotargs{:});
      annotaxis(fdepk.lat+0.005,'FDEP "K": met. \rightarrow',annotargs{:});
    end;
    if ( exist('ftpf1','var') )
      plot(ftpf1.lon,ftpf1.lat,plotargs{:});
      annotaxis(ftpf1.lat,'41114: Hs \rightarrow',annotargs{:});
    end;

    if ( exist('cnnf1','var') )
      plot(cnnf1.lon,cnnf1.lat,plotargs{:});
      annotaxis(cnnf1.lat,'41113: Hs \rightarrow',annotargs{:});
    end;
    if ( exist('canf1','var') )
      plot(canf1.lon,canf1.lat,plotargs{:});
      annotaxis(canf1.lat,'41009: met. \rightarrow',annotargs{:});
    end;

    if ( exist('sefcri','var') )
      plot(sefcri.dc1.lon,sefcri.dc1.lat,plotargs{:});
      plot(sefcri.dc3.lon,sefcri.dc3.lat,plotargs{:});
      annotaxis(sefcri.dc3.lat,'SEFCRI DC1/3 \rightarrow',annotargs{:});
      plot(sefcri.bc1.lon,sefcri.bc1.lat,plotargs{:});
      plot(sefcri.bc3.lon,sefcri.bc3.lat,plotargs{:});
      if ( blowupmap )
        annotaxis(sefcri.bc1.lat,'SEFCRI BC1 \rightarrow',annotargs{:});
        annotaxis(sefcri.bc3.lat,'SEFCRI BC3 \rightarrow',annotargs{:});
      else
        %annotaxis(sefcri.bc3.lat,'SEFCRI BC1/3 \rightarrow',annotargs{:});
        annotaxis(sefcri.bc3.lat+0.015,'SEFCRI BC1/3 \rightarrow',annotargs{:});
      end;

      plot(sefcri.pb1.lon,sefcri.pb1.lat,plotargs{:});
      plot(sefcri.pb2.lon,sefcri.pb2.lat,plotargs{:});
      if ( blowupmap )
        annotaxis(sefcri.pb1.lat,'SEFCRI PB1 \rightarrow',annotargs{:});
        annotaxis(sefcri.pb2.lat,'SEFCRI PB2 \rightarrow',annotargs{:});
      else
        annotaxis(sefcri.pb2.lat,'SEFCRI PB1/2 \rightarrow',annotargs{:});
      end;

      plot(sefcri.mc2.lon,sefcri.mc2.lat,plotargs{:});
      annotaxis(sefcri.mc2.lat,'SEFCRI MC2 \rightarrow',annotargs{:});
%%%%      if ( ~exist('fhix','var') || fhix > 1 )
        plot(sefcri.pela.lon,sefcri.pela.lat,plotargs{:});
        plot(sefcri.evca.lon,sefcri.evca.lat,plotargs{:});
        plot(sefcri.slib.lon,sefcri.slib.lat,plotargs{:});
%%%%      end;
      plot(sefcri.updb.lon,sefcri.updb.lat,plotargs{:});
      annotaxis(sefcri.updb.lat,'SEFCRI UpDB \rightarrow',annotargs{:});
    end; %if ( exist('sefcri','var') )

    if ( exist('fwc','var') )
      plot(fwc.TUCA.lon,fwc.TUCA.lat,plotargs{:});
      annotaxis(fwc.TUCA.lat,'FWC TUCA \rightarrow',annotargs{:});

      plot(fwc.Noula.lon,fwc.Noula.lat,plotargs{:});
      if ( blowupmap )
        annotaxis(fwc.Noula.lat,'FWC Noula \rightarrow',annotargs{:});
      else
        annotaxis(fwc.Noula.lat+0.020,'FWC Noula \rightarrow',annotargs{:});
      end;

      plot(fwc.Alpha.lon,fwc.Alpha.lat,plotargs{:});
      annotaxis(fwc.Alpha.lat,'FWC Alpha/Rodeo \rightarrow',annotargs{:});
      plot(fwc.Rodeo.lon,fwc.Rodeo.lat,plotargs{:});
    end;

    if ( exist('face','var') )
      plot(face.shallow.lon,face.shallow.lat,plotargs{:});
%%%%      if ( ~exist('fhix','var') || fhix > 1 )
        plot(face.tcm1.lon,face.tcm1.lat,plotargs{:});
        plot(face.tcm2.lon,face.tcm2.lat,plotargs{:});
%%%%      end;
      plot(face.deep.lon,face.deep.lat,plotargs{:});

      if ( ~isfield(face,'AOML_HW') )
        annotaxis(face.deep.lat,'FACE HWD moorings: u,v \rightarrow',annotargs{:});
      else
        plot(face.AOML_HW.lon,face.AOML_HW.lat,plotargs{:});
        plot(face.HanS_HW.lon,face.HanS_HW.lat,plotargs{:});
        if ( blowupmap )
          annotaxis(face.AOML_HW.lat,'FACE/H&S HW: u,v \rightarrow',annotargs{:});
        else
          annotaxis(face.AOML_HW.lat-0.050,'FACE/H&S HW: u,v \rightarrow',annotargs{:});
        end;
        plot(face.AOML_BR.lon,face.AOML_BR.lat,plotargs{:});
        plot(face.HanS_BR.lon,face.HanS_BR.lat,plotargs{:});
        if ( blowupmap )
          annotaxis(face.AOML_BR.lat,'FACE/H&S BR: u,v \rightarrow',annotargs{:});
        else
          annotaxis(face.AOML_BR.lat+0.020,'FACE/H&S BR: u,v \rightarrow',annotargs{:});
        end;
      end;

      if ( isfield(face,'brwd20') )
        plot(face.brwd20.lon,face.brwd20.lat,plotargs{:});
        plot(face.brwd40.lon,face.brwd40.lat,plotargs{:});
        plot(face.brwd100.lon,face.brwd100.lat,plotargs{:});
        plot(face.boca20.lon,face.boca20.lat,plotargs{:});
        plot(face.boca40.lon,face.boca40.lat,plotargs{:});
      end;
    end;

    %{ OLDER ALTERNATIVE NAME FOR THIS STRUCT
    if ( exist('jack','var') )
      plot(jack.shallow.lon,jack.shallow.lat,plotargs{:});
%%%%      if ( ~exist('fhix','var') || fhix > 1 )
        plot(jack.tcm1.lon,jack.tcm1.lat,plotargs{:});
        plot(jack.tcm2.lon,jack.tcm2.lat,plotargs{:});
%%%%      end;
      plot(jack.deep.lon,jack.deep.lat,plotargs{:});
      annotaxis(jack.deep.lat,'FACE HWD moorings: u,v \rightarrow',annotargs{:});
    end;
    %}

    if ( exist('sfomc','var') )
      plot(sfomc.ne_buoy.lon,sfomc.ne_buoy.lat,plotargs{:});
%%%%      if ( ~exist('fhix','var') || fhix > 1 )
        plot(sfomc.c_buoy.lon,sfomc.c_buoy.lat,plotargs{:});
        plot(sfomc.sw_buoy.lon,sfomc.sw_buoy.lat,plotargs{:});
        plot(sfomc.se_buoy.lon,sfomc.se_buoy.lat,plotargs{:});
        plot(sfomc.pier_cc.lon,sfomc.pier_cc.lat,plotargs{:});
%%%%      end;
      plot(sfomc.nw_w_btm.lon,sfomc.nw_w_btm.lat,plotargs{:});
      if ( blowupmap )
        annotaxis(sfomc.nw_w_btm.lat,'NSU NW/W,C,NE/E: Hs,u,v \rightarrow',annotargs{:});
        annotaxis(sfomc.sw_buoy.lat,'NSU SW buoy: u,v \rightarrow',annotargs{:});
      else
        annotaxis(sfomc.nw_w_btm.lat,'NSU moorings: Hs,u,v \rightarrow',annotargs{:});
      end;
    end;

