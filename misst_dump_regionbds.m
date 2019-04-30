function [xlon,xlat,dx] = misst_dump_regionbds(ireg)

  areg={'asam';'freef';'ecarib';'gbr'};
  dx=360/4096;

  if ( nargin < 1 || ireg < 1 || ireg > numel(areg) ); areg, return; end;

  disp(areg{ireg});
  ireg1(1)=2018;ireg1(2)=2986;ireg1(3)=3190;ireg1(4)=1529;
  jreg1(1)=687;jreg1(2)=1109;jreg1(3)=1029;jreg1(4)=653;
  ireg1=ireg1+8;jreg1=jreg1+8;
  [ireg1(ireg),jreg1(ireg)],
  xlon=([ireg1(ireg) ireg1(ireg)+256]-1)*dx+dx/2;
  xlat=([jreg1(ireg) jreg1(ireg)+256]-1)*dx-(90-dx/2);

  -(360 - xlon - (dx/2)),
  (xlat - (dx/2)),
  dx,

return;
