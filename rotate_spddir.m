function [newdr,u,v] = rotate_spddir(sp,dr,th)
%function [newdr,u,v] = rotate_spddir(sp,dr,th)
%
% Return a new vector NEWDR of directions in degrees, that is the result of
% rotating the vectors calculated by vectors SP and DR by TH degrees.  Calls
% SPDDIR_TO_UV, REORIENT_VECTORS, and UV_TO_DIR (v.), in that order.
%
% Last Saved Time-stamp: <Wed 2012-10-17 16:24:01 Eastern Daylight Time lew.gramer>

  [u,v]=spddir_to_uv(sp,dr);
  [u,v]=reorient_vectors(th,u,v);
  newdr=uv_to_dir(u,v);
  if ( nargout<2 )
    u=[];
    v=[];
    clear u v;
  end;

return;
