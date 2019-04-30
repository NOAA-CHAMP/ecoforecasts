1;

error('This is an old script');

stnam = 'fwyf1';

station = [];

if ( ~exist('datapath', 'var') || isempty(datapath) )
  datapath = './data';
end;

yrs = 1984:2008;
for ix = 1:length(yrs)
  yr = yrs(ix);
  fname = sprintf('%s/%sh%d.txt', datapath, lower(stnam), yr);
  if ( exist(fname, 'file') )
    station = load_ndbc_data(station, fname);
  end;
end;

clear datapath yrs ix yr fname;
