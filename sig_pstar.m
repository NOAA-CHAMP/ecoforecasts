function strs = sig_pstar(p_or_struct,more_stars)
%function strs = sig_pstar(p_or_struct,more_stars)
%
% Return a CHAR of zero or more asterisks corresponding to the significance
% of the numeric p-statistic P (or STRUCT with field .p). P <= 0.10: '.',
% 0.05: '*', 0.01: '**', 0.001: '***', 0.0001: '****'. If second arg is True
% (DEFAULT: False), continues appending up to '*******' for p<=1e-[5,6,7].
%
% Last Saved Time-stamp: <Sun 2019-02-24 14:43:42 Eastern Standard Time gramer>

  if ( isnumeric(p_or_struct) )
    p = p_or_struct;
  elseif ( isstruct(p_or_struct) )
    p = p_or_struct.p;
  end;
  if ( ~exist('more_stars') )
    more_stars = false;
  end;

  strs = '';
  if ( 0.05<p && p<=0.10 )
    strs='.';
  else
    if (p<=0.05);   strs=[strs '*']; end;
    if (p<=0.01);   strs=[strs '*']; end;
    if (p<=0.001);  strs=[strs '*']; end;
    if (p<=0.0001); strs=[strs '*']; end;
    if ( more_stars )
      if (p<=1e-5); strs=[strs '*']; end;
      if (p<=1e-6); strs=[strs '*']; end;
      if (p<=1e-7); strs=[strs '*']; end;
    end;
  end;

return;
