1;

x = importdata('data/mlrf1-1992-2002.xls');

dts = datenum(x.textdata(2:end,1));

t.date=dts; t.data = x.data(:,13); t.date(t.data<0)=[]; t.data(t.data<0)=[]; t.date(t.data>40)=[]; t.data(t.data>40)=[];
n.date=datenum(x.textdata(2:end,1)); n.data = x.data(:,14);
s.date=datenum(x.textdata(2:end,1)); s.data = x.data(:,12); s.date(s.data<0)=[]; s.data(s.data<0)=[];

scatter_fit_ts(t,mlrf1.ndbc_sea_t);
