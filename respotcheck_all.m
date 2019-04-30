1;

cfnms = {'30m','30m_NANMEAN_9_9_10','92m','92m_NANMEAN_3_3_4'};
cstnms={'CCAF1','CNDF1','CGRF1','CPIF1'};

for cfix = 1:numel(cfnms)
  fnm = cfnms{cfix};
  load(['FRT_depth_and_beta_',fnm,'_hc_range_depth.mat']);
  %for cstix=1:numel(cstnms)
  for cstix=3:3
    stnm=cstnms{cstix};
    fhix = ((cfix-1)*numel(cstnms))+1;
    fh = fmg(fhix);
    spotcheck_hc_range_field;
    appendtitlename([' ',textize(fnm)]);
    if ( cfix==1 ); ax{cstix}=axis; else; axis(ax{cstix}); end;
  end;  
  clear ans bet cs DBG dx flatix fnm h hch hcrng lat lon stnm
end;

clear cfnms cstnms cfix
