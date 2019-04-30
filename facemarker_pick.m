function mrkr = facemarker_pick(ix,mrkrs,nclrs)
%function mrkr = facemarker_pick(ix,mrkrs,nclrs)
%
% Return the (IX/NCLRS)th face Marker in the list (CHAR or cell array) of
% markers MRKRS (DEFAULT: {'-',':','-.','--'}). If IX > NUMEL(MRKRS), cycle
% through the markers repeatedly. See PLOT for the list of face markers. If
% IX is a CHAR string 'N', then the number of elements in MRKRS is returned.
%
% NOTE: The purpose of this function is to be used in concert with COLOR_PICK
% (v.) to ensure a unique color/marker combo for up to NCLRS*NUMEL(MRKRS)
% different indices IX, as in the following example:
%
%  >> for ix=1:24; plot(1:3,rand([1,3]),'Color',color_pick(ix),'Marker',facemarker_pick(ix)); end;
%
% Last Saved Time-stamp: <Mon 2017-05-01 15:23:51 Eastern Daylight Time gramer>

  if ( ~exist('mrkrs','var') || isempty(mrkrs) )
    mrkrs = {'.','o','x','+','*','s','d','v','^','<','>','p','h'};
  end;
  if ( ~exist('nclrs','var') || isempty(nclrs) )
    % This DEFAULT number should match the length of the DEFAULT color array in COLOR_PICK (v.)
    %nclrs = 6;
    nclrs = color_pick('ncolors');
  end;
  if ( ~iscell(mrkrs) )
    if ( ischar(mrkrs) )
      mrkrs = cellstr(mrkrs(:));
    else
      error('Arg MRKRS must be a char or cell array of face Marker strings');
    end;
  end;
  if ( ischar(ix) && strncmpi(ix,'N',1) )
    mrkr = numel(mrkrs);
  else
    mrkr = mrkrs{mod(floor((ix-1)/nclrs),numel(mrkrs))+1};
  end;

return;
