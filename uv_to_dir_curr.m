function cdir = uv_to_dir_curr(u, v)
%function cdir = uv_to_dir_curr(u, v)
%
% Convert u and v vector components into a direction in degrees True.
% NOTE: This version of 'uv-to-dir' is coded for OCEAN CURRENTS: "dir" here
% means "target direction", NOT "source direction" as it would for winds.
% SEE ALSO: UV_TO_DIR
%
% Last Saved Time-stamp: <Sun 2010-07-25 14:28:01 Eastern Daylight Time gramer>

    cdir = uv_to_dir(u, v);

    cdir = cdir - 180;
    negix = find(cdir < 0);
    cdir(negix) = 360 + cdir(negix);

return;
