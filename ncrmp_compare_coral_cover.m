1;

%A = ncrmp.FGBNMS_2013;		B = ncrmp.FGBNMS_2015;
%A = ncrmp.SEFCRI_2016;		B = ncrmp.SEFCRI_2014;
%A = ncrmp.FLK_2014;		B = ncrmp.FLK_2016;
%A = ncrmp.TortugasMarq_2014;	B = ncrmp.TortugasMarq_2016;
A = ncrmp.PRICO_2014;		B = ncrmp.PRICO_2016;
%A = ncrmp.USVI_2015;		B = ncrmp.USVI_2017;
%%A = ncrmp.USVI_2015;		B = ncrmp.USVI_2013; % Just northern islands!

clear mindx minix
Bix = 1:numel(B.ulon);
for ix=1:numel(B.ulon)
  [mindx(ix),minix(ix)] = min(distance_wgs84(B.ulat(ix),B.ulon(ix),A.ulat,A.ulon));
end;

% Filter out widely separated pairs
badix = find(mindx>1);
mindx(badix) = [];
minix(badix) = [];
Bix(badix) = [];

% fmg; plot(B.ulon(Bix),B.ulat(Bix),'r.',A.ulon(minix),A.ulat(minix),'b.'); daspect([1,cosd(B.ulat(Bix(1))),1]); titlename(textize([A.basename,' vs. ',B.basename]));

scatter_fit(A.utc(minix),B.utc(Bix),textize(A.basename),textize(B.basename),[],[],[],true,{'Location','Best'}); axis([0,100,0,100]);

%scatter_fit(mindx,abs(A.utc(minix)-B.utc(Bix)),textize(['DX(',A.basename,',',B.basename,')']),textize(['abs(',A.basename,'-',B.basename,')'])); axis([0,100,0,100]);


clear neardx nearix
all_ixen = 1:numel(A.ulon);
for ix=1:numel(A.ulon)
  ixen = all_ixen(all_ixen ~= ix);
  [neardx(ix),nearix(ix)] = min(distance_wgs84(A.ulat(ix),A.ulon(ix),A.ulat(ixen),A.ulon(ixen)));
end;
clear all_ixen ixen

fmg; hist(neardx,100); title('Min DX');
scatter_fit(neardx,abs(A.utc-A.utc(nearix)),'Nearest DX','Total Cvr % diff.');

%clear ans Bix badix ix mindx minix neardx nearix
