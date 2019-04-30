1;

%doPrint = false;
doPrint = true;

% Linear R2 ~ 0.58, Exp(-1/3) R2 ~ 0.92! (N=30)
doAccDist = false; HsRange=[]; fillBathGaps = false; plot_range_vs_control
 bat=[]; bath=[]; kwbat=[]; kwbath=[]; stns=[]; clear -regexp ^[a-cegi-z]*; pack

% For 7-31 to 10-1, linear R2 ~ 0.74! (N=44)
doAccDist = true; HsRange=[]; fillBathGaps = false; ctlvar = 'raw_nwps_sigwavehgt'; plot_range_vs_control
 bat=[]; bath=[]; kwbat=[]; kwbath=[]; stns=[]; clear -regexp ^[a-cegi-z]*; pack

% For 7-31 to 10-1, linear R2 ~ 0.56, exp(-1/3) R2 ~ 0.91! (N=28)
HsRange=[0,0.5]; fillBathGaps = false; plot_range_vs_control
 bat=[]; bath=[]; kwbat=[]; kwbath=[]; stns=[]; clear -regexp ^[a-cegi-z]*; pack


doAccDist = true; HsRange=[0,0.5]; fillBathGaps = false; plot_range_vs_control

persub = @ts_jas; HsRange=[]; fillBathGaps = false; ctlvar = 'raw_nwps_sigwavehgt'; plot_range_vs_control
 bat=[]; bath=[]; kwbat=[]; kwbath=[]; stns=[]; clear -regexp ^[a-cegi-oq-z]*; pack

persub = @ts_jas; HsRange=[0,0.4]; fillBathGaps = false; plot_range_vs_control
persub = @ts_jas; HsRange=[0,0.5]; fillBathGaps = false; plot_range_vs_control
 bat=[]; bath=[]; kwbat=[]; kwbath=[]; stns=[]; clear -regexp ^[a-cegi-oq-z]*; pack


%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doPrint = false;

% For 7-31 to 10-1, linear R2 ~ 0.49, exp(-1/3) R2 ~ 0.44! (N=32)
% (Setting fillBathGaps adds a "Channel" site! But also a "Shoals", and...?)
HsRange=[0,0.5]; fillBathGaps = true; plot_range_vs_control
 bat=[]; bath=[]; kwbat=[]; kwbath=[]; stns=[]; clear -regexp ^[a-cegi-z]*; pack

% For 7-31 to 10-31, linear R2 ~ 0.77! (N=10)
HsRange=[0,0.45]; fillBathGaps = false; plot_range_vs_control
 bat=[]; bath=[]; kwbat=[]; kwbath=[]; stns=[]; clear -regexp ^[a-cegi-z]*; pack

% For 7-31 to 10-1, linear R2 ~ 0.79! (N=36)
HsRange=[0,0.5]; fillBathGaps = true; ctlvar = 'raw_nwps_sigwavehgt'; plot_range_vs_control
 bat=[]; bath=[]; kwbat=[]; kwbath=[]; stns=[]; clear -regexp ^[a-cegi-z]*; pack

% For 7-31 to 10-1, linear R2 ~ 0.74! (N=44)
HsRange=[]; fillBathGaps = true; ctlvar = 'raw_nwps_sigwavehgt'; plot_range_vs_control
 bat=[]; bath=[]; kwbat=[]; kwbath=[]; stns=[]; clear -regexp ^[a-cegi-z]*; pack


%% Experiment: Look at HIGH wave-energy sites instead...
HsRange = [0.3,+inf];
%fillBathGaps = true;
fillBathGaps = false;
 ctlvar = 'ngdc_depth'; plot_range_vs_control
 bat=[]; bath=[]; kwbat=[]; kwbath=[]; stns=[]; clear -regexp ^[a-cegi-z]*; pack
 ctlvar = {'ngdc_depth',@nanmean,4,4,5}; plot_range_vs_control
 bat=[]; bath=[]; kwbat=[]; kwbath=[]; stns=[]; clear -regexp ^[a-cegi-z]*; pack
 ctlvar = 'ngdc_beta'; plot_range_vs_control
 bat=[]; bath=[]; kwbat=[]; kwbath=[]; stns=[]; clear -regexp ^[a-cegi-z]*; pack
 ctlvar = 'raw_nwps_sigswellhgt'; plot_range_vs_control
 bat=[]; bath=[]; kwbat=[]; kwbath=[]; stns=[]; clear -regexp ^[a-cegi-z]*; pack
 % ctlvar = {'ngdc_beta',@nanmean,3,3,5}; plot_range_vs_control
 % bat=[]; bath=[]; kwbat=[]; kwbath=[]; stns=[]; clear -regexp ^[a-cegi-z]*; pack
 ctlvar = {'ngdc_beta',@nanmean,4,4,5}; plot_range_vs_control
 bat=[]; bath=[]; kwbat=[]; kwbath=[]; stns=[]; clear -regexp ^[a-cegi-z]*; pack
 % ctlvar = {'ngdc_beta',@nanmean,5,5,5}; plot_range_vs_control
 % bat=[]; bath=[]; kwbat=[]; kwbath=[]; stns=[]; clear -regexp ^[a-cegi-z]*; pack

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
