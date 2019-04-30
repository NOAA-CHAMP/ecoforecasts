function tOrF = isOctave
%function tOrF = isOctave
% Return TRUE if we are running in Octave (vs. MATLAB): Lew.Gramer@noaa.gov

  persistent isOct;

  if ( isempty(isOct) )
    v = ver;
    if ( strcmpi(v(1).Name,'Octave') )
      isOct = true;
    else
      isOct = false;
    end;
  end;

  tOrF = isOct;

return;
