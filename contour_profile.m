function contour_profile(pstr,fld,ix,varargin)
%function contour_profile(pstr,fld,ix,varargin)
%
% Use CONTOURF to plot profile FLD using dates PSTR.date and depths
% in PSTR.(DEPFLD) (DEFAULT: PSTR.depths).  DEFAULT FLD is STRUCT field
% FLDSTR.field. If IX is given, plot SQUEEZE(FLD(ix,:,:)), unless IX is a
% FUNCTION_HANDLE (v.): in that case, plot SQUEEZE(IX(FLD)). All other
% arguments are passed through to CONTOURF as its args 4-N.
%
% Last Saved Time-stamp: <Fri 2016-03-18 14:49:12 Eastern Daylight Time gramer>

  if ( ~exist('fld','var') || isempty(fld) )
    fld = pstr.prof;
  elseif ( ischar(fld) )
    fld = pstr.(fld);
  end;
  if ( exist('ix','var') && ~isempty(ix) )
    if ( isa(ix,'function_handle') )
      fld = ix(fld);
    else
      fld = fld(ix,:,:);
    end;
  end;
  if ( ~exist('dfld','var') || isempty(dfld) )
    dfld = 'depths';
  end;

  fld = squeeze(fld);

  contourf(pstr.date,pstr.(dfld),fld',varargin{:});
  if ( numel(varargin) > 0 && numel(varargin{1}) > 2 )
    set(gca,'CLim',[varargin{1}(1),varargin{1}(end)]);
  end;
  colorbar;
  if ( exist('datetick3') )
    disp('datetick3');
    datetick3;
  elseif ( exist('datetick2_5') )
    datetick2_5;
  elseif ( exist('datetick2') )
    datetick2;
  else
    datetick;
  end;

return;
