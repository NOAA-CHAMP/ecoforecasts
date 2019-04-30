function nms = get_monthnames
%function nms = get_monthnames
% Return a char array of the 12 month names for the current locale
% Last Saved Time-stamp: <Fri 2014-10-24 13:18:54 Eastern Daylight Time gramer>
  nms = datestr(datenum(1,1:12,1),'mmm');
return;
