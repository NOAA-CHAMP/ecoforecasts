1;


for qix=1:4;
  % Net sea-surface heat flux, Q_0
  q0 = q0s(qix);
  if ( strcmp(subrgn,'UK') )
    disp('Hit Enter to zoom to sub-region'); pause;
    figure(qix);
    axis([-80.21,-80.135,25.33,25.40]); daspect([1,cosd(25),1]); view(2);
    colorbar('off'); colorbar('Location','EastOutside');
    if ( doPrint )
      print('-dpng',fullfile(figspath,sprintf('%s_%s_extreme_%dW_BLOWUP.png',figbasename,subrgn,abs(q0s(qix)))));
    end;
  end;
end;
