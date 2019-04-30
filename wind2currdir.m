function wdir = wind2currdir(udir)

  wdir = udir + 180;
  wdir(wdir >= 360.0) = wdir(wdir >= 360.0) - 360;

return;
