function sst = read_misst(fname,lonlen,latlen)
%function sst = read_misst(fname,lonlen,latlen)
%
% Read a raw MISST binary file named FNAME, laid out with LONLEN zonal and
% LATLEN merid. gridpoints. Return LATLENxLONLEN matrix of SST [DOUBLE, oC].
%
%% Two sample calls to READ_MISST that should produce identical figures
% sst = read_misst('data/misst/mw_ir.freef.fusion.2009.365.v02',256,256);
% figure; maxigraph; contourf(sst,[15:28]); colorbar; title('FReef');
%
% ox=2986; oy=1109; ox=ox+8+1; oy=oy+8+1;
% gsst = read_misst('data/misst/mw_ir.fusion.2009.365.v02',4096,2048);
% figure; maxigraph; contourf(gsst(oy:oy+255,ox:ox+255),[15:28]); colorbar; title('Global');
%
% Last Saved Time-stamp: <Wed 2011-02-09 13:36:06  lew.gramer>

  fid = fopen(fname, 'r');
  if ( fid < 0 || feof(fid) || ~isempty(ferror(fid)) )
    error('Unable to open "%s"!',fname);
  end;
  % for ix=1:lonlen
  %   sst(ix,1:latlen) = cast(fread(fid,latlen,'uchar'),'double');
  % end;
  dat = cast(fread(fid,(lonlen*latlen),'uchar'),'double');
  fclose(fid);

  sst = reshape(dat', [lonlen latlen])';
  dat = []; clear dat;

  %% 250 is also a mask value apparently!
  % sst(sst > 250) = nan;
  sst(sst >= 250) = nan;
  sst = (sst .* 0.15) - 3;

return;
