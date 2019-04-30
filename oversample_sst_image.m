function [osst,omask,sst,mask] = oversample_sst_image(fpath,subs,isAnom,gapMethod,interpMethod)
%function [osst,omask,sst,mask] = oversample_sst_image(fpath,subs,isAnom,gapMethod,interpMethod)
%
% Return Sea Surface Temperature field from file FPATH, oversampled to
% fractional resolution SUBS using NaN filling and INTERP2 (v.). For example,
% OVERSAMPLE_SST_IMAGE(F,0.2) on a 4km SST image in F gives 800m resolution.
% Arg GAPMETHOD controls gap filling on the source image: it may be 'raw' or
% 'none' to do no filling (DEFAULT), or NANMEAN, NANMEDIAN or other function
% to do 2x2 pixel averages for gap-filling. Arg INTERPMETHOD if given is
% passed in to INTERP2 (v.). See GET_SST_IMAGE for the use of arg ISANOM.
%
% Oversampled SST field and oversampled land-mask are both returned, along
% with original-resolution SST field and land-mask.
%
% Last Saved Time-stamp: <Mon 2011-03-07 19:02:09  lew.gramer>

  if ( ~exist('isAnom','var') || isempty(isAnom) )
    isAnom = false;
  end;
  if ( ~exist('gapMethod','var') || isempty(gapMethod) )
    gapMethod = 'none';
  end;
  if ( ~exist('interpMethod','var') || isempty(interpMethod) )
    interpMethod = 'linear';
  end;

  % Assume interpolation grid is always equal-spaced and monotonic
  if ( interpMethod(1) ~= '*' )
    interpMethod = ['*' interpMethod];
  end;


  [sst,mask] = read_sst_image(fpath,isAnom);

  psst = sst;
  switch ( gapMethod )
   case {'','none','raw'},
   otherwise,
    [nanix,nanjx]=find(~isfinite(psst));
    bdrix=find(nanix==1|nanix==size(psst,1)|nanjx==1|nanjx==size(psst,2));
    nanix(bdrix)=[]; nanjx(bdrix)=[];
    for ix=1:length(nanix)
      mysst=sst(nanix(ix)-1:nanix(ix)+1,nanjx(ix)-1:nanjx(ix)+1);
      psst(nanix(ix),nanjx(ix)) = feval(gapMethod,mysst(:));
    end;
  end;

  osst = interp2(psst,[1:subs:size(psst,2)]',[1:subs:size(psst,1)],interpMethod);

  % Overlay land mask on oversampled grid in a reasonably conservative way
  omask = logical(interp2(mask,[1:subs:size(mask,2)]',[1:subs:size(mask,1)],'*linear'));
  osst(omask) = nan;

  if ( nargout<2 )
    omask=[]; clear omask;
  elseif ( nargout<3 )
    sst=[]; clear sst;
  elseif ( nargout<4 )
    mask=[]; clear mask;
  end;

return;
