function [rgs,ixes,allixes] = find_date_ranges(dts,maxgap)
%function [rgs,ixes,allixes] = find_date_ranges(dts,maxgap)
%
% Find beginning and ending of all "contiguous" sequences in DTS, a vector of
% DATENUMs. "Contiguous" sequences are those having no gap longer than MAXGAP
% (DEFAULT: +Inf). If no return value, pretty-print contiguous date ranges.
% Otherwise, return RGS a 2xN matrix of DATENUMs for contiguous ranges, IXES
% a 2xN matrix of gap indices, and ALLIXES the vector of all good indices.
%
% Last Saved Time-stamp: <Fri 2017-05-26 14:23:14 Eastern Daylight Time gramer>

  if ( ~exist('maxgap','var') || isempty(maxgap) )
    maxgap = +Inf;
  end;

  dtsdif = diff(dts(:));
  gapix = [1,find(0>dtsdif|dtsdif>maxgap)'+1,length(dts)+1];

  allixes = [];
  for ix = 1:length(gapix)-1
    allixes = [allixes,gapix(ix):gapix(ix+1)-1];
    ixes(1:2,ix) = [gapix(ix),gapix(ix+1)-1];
    rgs(1:2,ix) = dts(ixes(1:2,ix));
  end;

  if ( nargout < 1 )
    totalgaps = 0;
    totaldays = 0;
    for ix = 1:length(gapix)-1
      difdt = rgs(2,ix) - rgs(1,ix);
      numdy = floor(difdt); numhr = 24*(difdt-floor(difdt));
      numdys = sprintf(' [% 4dd% 3dh]',floor(numdy),round(numhr));
      if ( ix < size(rgs,2) )
        gapdy = rgs(1,ix+1)-rgs(2,ix);
        totalgaps = totalgaps + gapdy;
        gaphr = 24*(gapdy-floor(gapdy));
        gapdys = sprintf(' (% 4dd% 3dh)',floor(gapdy),round(gaphr));
      else
        gapdys = '';
      end;
      disp([ sprintf('% 8d',ixes(1,ix)) ' ' datestr(rgs(1,ix),0) ' - ' ...
             sprintf('% 8d',ixes(2,ix)) ' ' datestr(rgs(2,ix),0) numdys gapdys]);

      daydy = rgs(2,ix)-rgs(1,ix);
      totaldays = totaldays + daydy;
    end;

    if ( totalgaps > 0 )
      gaphr = 24*(totalgaps-floor(totalgaps));
      gappct = 100*totalgaps/(dts(end)-dts(1));
      gapdys = sprintf(' % 5dd% 3dh  (% 4.2f%%)',floor(totalgaps),round(gaphr),gappct);
      disp(['Total gaps: ',gapdys]);
    end;

    dayhr = 24*(totaldays-floor(totaldays));
    daypct = 100*totaldays/(dts(end)-dts(1));
    daydys = sprintf(' % 5dd% 3dh  (% 4.2f%%)',floor(totaldays),round(dayhr),daypct);
    disp(['Total days: ',daydys]);

    allixes=[]; ixes=[]; rgs=[]; clear allixes ixes rgs;
  end;

return;
