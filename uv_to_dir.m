function wdir = uv_to_dir(u, v)
%function wdir = uv_to_dir(u, v)
%
% Convert u and v vector components into a direction in degrees True
% NOTE: This version of 'uv-to-dir' is coded for WINDS: "dir" here means
% "source direction", NOT "target direction" as it would for ocean currents.
% SEE ALSO: UV_TO_DIR_CURR
%
% Last Saved Time-stamp: <Tue 2015-01-06 23:05:22 Eastern Standard Time gramer>

    wdir = repmat(nan, size(u));

    spd = uv_to_spd(u, v);

    % Handle direction for "zero wind" specially
    bpidx = find((u>0 & v>=0) | (u>=0 & v>0));
    vpidx = find((u<0 & v>=0) | (u<=0 & v>0));
    upidx = find((u>0 & v<=0) | (u>=0 & v<0));
    npidx = find((u<0 & v<=0) | (u<=0 & v<0));
    b0idx = find(u==0 & v==0);

    wdir(bpidx) = ( 180 + asind(u(bpidx) ./ spd(bpidx)) );
    wdir(vpidx) = ( 180 + asind(u(vpidx) ./ spd(vpidx)) );
    wdir(upidx) = ( 360 - asind(u(upidx) ./ spd(upidx)) );
    wdir(npidx) = ( - asind(u(npidx) ./ spd(npidx)));
    % Make sure this one happens last
    wdir(b0idx) = 0;

return;
