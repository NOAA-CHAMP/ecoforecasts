function [spd,dir] = uv_to_spddir_curr(u, v)
%function [spd,dir] = uv_to_spddir_curr(u, v)
% Convert u and v vector components into a magnitude (a speed) and TO direction
% Last Saved Time-stamp: <Sun 2017-11-19 15:39:15 Eastern Standard Time gramer>

  spd = uv_to_spd(u,v);
  dir = uv_to_dir_curr(u,v);

return;
