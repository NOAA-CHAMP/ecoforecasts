function strs = get_season_name(seas,longform)
%function strs = get_season_name(seas,longform)
% Get "season name" for each season number (1-4) SEAS. "Season number" is the
% equivalent of a month triad: JFM, AMJ, JAS, OND (integers 1-4, resp.) If
% LONGFORM is True, return e.g., 'Jan-Mar' in place of 'JFM' for US-English.
% Note LONGFORM may also be a LOGICAL array with NUMEL == NUMEL(SEAS).
% 
% Last Saved Time-stamp: <Wed 2018-08-29 16:56:53 Eastern Daylight Time gramer>

  strs = repmat(' ',size(seas));
  if ( ~exist('longform','var') || isempty(longform) )
    longform = false;
  end;

  % Make strings Locale-specific
  %season_names      = ['JFM    ';'AMJ    ';'JAS    ';'OND    '];
  %long_season_names = ['Jan-Mar';'Apr-Jun';'Jul-Sep';'Oct-Dec'];
  monn = datenum(0,1:12,1);
  mons = datestr(monn,'mmm');
  season_names = [mons(1:3,1)','    ';mons(4:6,1)','    ';mons(7:9,1)','    ';mons(10:12,1)','    '];
  long_season_names = [mons(1,:),'-',mons(3,:);mons(4,:),'-',mons(6,:);mons(7,:),'-',mons(9,:);mons(10,:),'-',mons(12,:)];

  if ( isscalar(longform) )
    if ( longform )
      strs = long_season_names(seas,:);
    else
      strs = season_names(seas,:);
    end;
  elseif ( numel(longform) == numel(seas) )
    strsz = size(long_season_names,2);
    strs(find(longform),1:strsz) = char(long_season_names(seas(find(longform)),:));
    strs(find(~longform),1:strsz) = char(season_names(seas(find(~longform)),:));
  else
    error('LONGFORM must be scalar or have NUMEL == NUMEL(SEAS)');
  end;

return;
