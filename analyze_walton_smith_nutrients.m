1;


%x = importdata('../../coral/SFP/walton_smith_database.xls');
x = importdata(get_ecoforecasts_path('data/walton_smith_database.xls'));

{x.textdata.Sheet1{1,:}}

goodix = find(~cellfun(@isempty,x.textdata.Sheet1(:,5)));
goodix(1) = [];
dts = datenum( {x.textdata.Sheet1{goodix,5}} );

z = x.data.Sheet1(goodix-1,4);
T1 = x.data.Sheet1(goodix-1,7);
N = x.data.Sheet1(goodix-1,18);
Si = x.data.Sheet1(goodix-1,19);
P = x.data.Sheet1(goodix-1,20);

scatter_fit(T1(T1<23&N>1),N(T1<23&N>1),'T','NO3+NO2');
titlename('Walton Smith, all years - AOML South Florida Program');
% print('-dtiff','figs/walton-smith-T-vs-N.tif');
% print('-dpng','figs/walton-smith-T-vs-N.png');
% print('-dtiff','../../coral/CRCP/Upwelling/CoRIS/walton-smith-T-vs-N.tif');
% print('-dpng','../../coral/CRCP/Upwelling/CoRIS/walton-smith-T-vs-N.png');


scatter_fit(T1(T1<23&P>0.05),P(T1<23&P>0.05),'T','PO4');
titlename('Walton Smith, all years - AOML South Florida Program');
% print('-dtiff','figs/walton-smith-T-vs-P.tif');
% print('-dpng','figs/walton-smith-T-vs-P.png');
% print('-dtiff','../../coral/CRCP/Upwelling/CoRIS/walton-smith-T-vs-P.tif');
% print('-dpng','../../coral/CRCP/Upwelling/CoRIS/walton-smith-T-vs-P.png');


%zcutoff=50; scatter_fit(T1(z>zcutoff),-z(z>zcutoff),'T',sprintf('z(z<-%g)',zcutoff)), ylim([-600,0])

%dtcrit=(get_year(dts)>=1997 & (get_season(dts)==2|get_season(dts)==3)); dtcritstr = 'Spring-Summer';
%dtcrit=(get_year(dts)>=2000 & ismember(get_month(dts),[4:11])); dtcritstr = 'Apr-Nov';
%dtcrit=(get_year(dts)>=1997 & ismember(get_month(dts),[4:11])); dtcritstr = 'Apr-Nov';
%dtcrit=(get_year(dts)>=1997 & ismember(get_month(dts),[4:10])); dtcritstr = 'Apr-Oct';

ztop = 10;
%ztop = 0;
Tbtm = 12;
Ttop = 24;

dtcrit=(get_year(dts)>=1997 & (get_season(dts)==2|get_season(dts)==3)); dtcritstr = 'Spring-Summer';

Nbtm = 1;
scatter_fit(T1(Tbtm<=T1&T1<=Ttop&z>ztop&N>=Nbtm&dtcrit),...
            N(Tbtm<=T1&T1<=Ttop&z>ztop&N>=Nbtm&dtcrit),...
            sprintf('%s T(%g<=T<=%g)',dtcritstr,Tbtm,Ttop),...
            sprintf('N(N>=%g,z<-%g)',Nbtm,ztop)),
ylim([0,26]);
legh=legend('Location','NorthEast');
%print('-dpng',sprintf('../../coral/CRCP/Upwelling/CoRIS/walton-smith-T-vs-N-%g-%g-%g-%g-%s.png',...
%                      ztop,Tbtm,Ttop,Nbtm,dtcritstr));
titlename('Sea temperature vs. [NO_3^-+NO_2^-]: south and west FRT');
set(gca,'LineWidth',1.5); set(get(gca,'Child'),'LineWidth',2); set(get(gca,'Child'),'MarkerSize',12);
set(legh,'FontSize',20);
print('-dpng',fullfile(get_coral_path,'CRCP','Upwelling','CoRIS',...
                       sprintf('walton-smith-T-vs-N-%g-%g-%g-%g-%s.png',...
                               ztop,Tbtm,Ttop,Nbtm,dtcritstr)));


%dtcrit=(get_year(dts)>=1997 & ismember(get_month(dts),[4:10])); dtcritstr = 'Apr-Oct';

%Pbtm = 0;
Pbtm = 0.0625;
scatter_fit(T1(Tbtm<=T1&T1<=Ttop&z>ztop&P>=Pbtm&dtcrit),...
            P(Tbtm<=T1&T1<=Ttop&z>ztop&P>=Pbtm&dtcrit),...
            sprintf('%s T(%g<=T<=%g)',dtcritstr,Tbtm,Ttop),...
            sprintf('P(P>=%g,z<-%g)',Pbtm,ztop)),
ylim([0,2]);
legh=legend('Location','NorthEast');
% print('-dpng',sprintf('../../coral/CRCP/Upwelling/CoRIS/walton-smith-T-vs-P-%g-%g-%g-%g-%s.png',...
%                       ztop,Tbtm,Ttop,Pbtm,dtcritstr));
titlename('Sea temperature vs. [P]: south and west FRT');
set(gca,'LineWidth',1.5); set(get(gca,'Child'),'LineWidth',2); set(get(gca,'Child'),'MarkerSize',12);
set(legh,'FontSize',20);
print('-dpng',fullfile(get_coral_path,'CRCP','Upwelling','CoRIS',...
                       sprintf('walton-smith-T-vs-P-%g-%g-%g-%g-%s.png',...
                               ztop,Tbtm,Ttop,Pbtm,dtcritstr)));


%Sibtm = 0;
Sibtm = 0.0625;
scatter_fit(T1(Tbtm<=T1&T1<=Ttop&z>ztop&Si>=Sibtm&dtcrit),...
            Si(Tbtm<=T1&T1<=Ttop&z>ztop&Si>=Sibtm&dtcrit),...
            sprintf('%s T(%g<=T<=%g)',dtcritstr,Tbtm,Ttop),...
            sprintf('Si(Si>=%g,z<-%g)',Sibtm,ztop)),
%ylim([0,2]);
ylim([0,10]);
legh=legend('Location','NorthEast');
% print('-dpng',sprintf('../../coral/CRCP/Upwelling/CoRIS/walton-smith-T-vs-Si-%g-%g-%g-%g-%s.png',...
%                       ztop,Tbtm,Ttop,Sibtm,dtcritstr));
titlename('Sea temperature vs. [Si]: south and west FRT');
set(gca,'LineWidth',1.5); set(get(gca,'Child'),'LineWidth',2); set(get(gca,'Child'),'MarkerSize',12);
set(legh,'FontSize',20);
print('-dpng',fullfile(get_coral_path,'CRCP','Upwelling','CoRIS',...
                       sprintf('walton-smith-T-vs-Si-%g-%g-%g-%g-%s.png',...
                               ztop,Tbtm,Ttop,Sibtm,dtcritstr)));

%clear x goodix dts z T1 N Si P ztop Tbtm Ttop dtcrit dtcritstr Nbtm Pbtm Sibtm


% Where is the mixed layer?


