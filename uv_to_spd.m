function spd = uv_to_spd(u, v)
%function spd = uv_to_spd(u, v)
% Convert u and v vector components into a magnitude (a speed)
% Last Saved Time-stamp: <Sun 2017-11-19 15:38:42 Eastern Standard Time gramer>

    spd = sqrt((u.*u) + (v.*v));

return;
