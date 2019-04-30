function styl = linestyle_pick(ix,styls,nclrs)
%function styl = linestyle_pick(ix,styls,nclrs)
%
% Return the (IX/NCLRS)th LineStyle in the list (CHAR or cell array) of
% styles STYLS (DEFAULT: {'-',':','-.','--'}). If IX > NUMEL(STYLS), cycle
% through the styles repeatedly. See PLOT for the list of styles.  If IX is a
% CHAR string 'N', then the number of elements in STYLS is returned.
%
% NOTE: The purpose of this function is to be used in concert with COLOR_PICK
% (v.) to ensure a unique color/style combo for up to NCLRS*NUMEL(STYLS)
% different indices IX, as in the following example:
%
%  >> for ix=1:24; plot(1:3,rand([1,3]),'Color',color_pick(ix),'LineStyle',linestyle_pick(ix)); end;
%
% Last Saved Time-stamp: <Mon 2017-05-01 15:25:19 Eastern Daylight Time gramer>

  if ( ~exist('styls','var') || isempty(styls) )
    styls = {'-',':','-.','--'};
  end;
  if ( ~exist('nclrs','var') || isempty(nclrs) )
    % This DEFAULT number should match the length of the DEFAULT color array in COLOR_PICK (v.)
    nclrs = 6;
  end;
  if ( ~iscell(styls) )
    if ( ischar(styls) )
      styls = cellstr(styls(:));
    else
      error('Arg STYLS must be a char or cell array of LineStyle strings');
    end;
  end;
  if ( ischar(ix) && strncmpi(ix,'N',1) )
    styl = numel(styls);
  else
    styl = styls{mod(floor((ix-1)/nclrs),numel(styls))+1};
  end;

return;
