1;

s = tdfread('1606021_SL  10_(0)_CR.txt',',');
dts = datenum(s.ISO_8601_Time,'yyyy-mm-ddTHH:MM:SS');
fmg; plot(dts,s.Velocity0x2DN_0x28cm0x2Fs0x29,dts,s.Velocity0x2DE_0x28cm0x2Fs0x29); datetick3; legend('S-N','W-E');
