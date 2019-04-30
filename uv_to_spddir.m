function [spd,dir] = uv_to_spddir(u, v)
%function [spd,dir] = uv_to_spddir(u, v)
% Convert u and v vector components into a magnitude (a speed) and direction
% Last Saved Time-stamp: <Sun 2017-11-19 15:38:58 Eastern Standard Time gramer>

  spd = uv_to_spd(u,v);
  dir = uv_to_dir(u,v);

return;
