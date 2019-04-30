1;
% Load daily Florida Current transport estimate (Sv) from 27N cable.
% To update cable transport file, run: d:/ecoforecasts/data/FC/get_FC_daily.csh
% 
%  "For further information [on Cable Transport] please contact
%  Dr. Christopher Meinen
%  Christopher.Meinen@noaa.gov
%  (305) 361-4355 "
%
% Lew.Gramer@noaa.gov

FC = load('d:/ecoforecasts/data/FC/FC_cable_transport.dat');
fc.date = datenum(FC(:,[1,2,3]));
fc.data = FC(:,4);
fc.qc = FC(:,5);
clear FC

qc_fc.date = fc.date;
qc_fc.data = fc.data;
qc_fc.data(fc.qc > 0) = nan;


