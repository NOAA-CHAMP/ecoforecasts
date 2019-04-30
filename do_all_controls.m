1;

doPrint = true;

%% SEAFLOOR DEPTH
persub = @ts_jfm; doWaves=false; ctlvar={'ngdc_depth',@nanmean,9,9,10}; doAccDist=false; plot_range_vs_control
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; diary off; more on; pack
persub = @ts_jas; doWaves=false; ctlvar={'ngdc_depth',@nanmean,9,9,10}; doAccDist=false; plot_range_vs_control
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; diary off; more on; pack
persub = @ts_jfm; doWaves=false; ctlvar={'ngdc_depth',@nanmean,9,9,10}; doAccDist=true; plot_range_vs_control
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; diary off; more on; pack
persub = @ts_jas; doWaves=false; ctlvar={'ngdc_depth',@nanmean,9,9,10}; doAccDist=true; plot_range_vs_control
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; diary off; more on; pack

persub = @ts_jfm; HsRange=[0,0.5]; ctlvar={'ngdc_depth',@nanmean,9,9,10}; doAccDist=false; plot_range_vs_control
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; diary off; more on; pack
persub = @ts_jas; HsRange=[0,0.5]; ctlvar={'ngdc_depth',@nanmean,9,9,10}; doAccDist=false; plot_range_vs_control
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; diary off; more on; pack
persub = @ts_jfm; HsRange=[0,0.5]; ctlvar={'ngdc_depth',@nanmean,9,9,10}; doAccDist=true; plot_range_vs_control
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; diary off; more on; pack
persub = @ts_jas; HsRange=[0,0.5]; ctlvar={'ngdc_depth',@nanmean,9,9,10}; doAccDist=true; plot_range_vs_control
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; diary off; more on; pack

%% WAVE EXPOSURE
persub = @ts_jfm; doWaves=false; ctlvar={'ngdc_depth',@nanmean,9,9,10}; doAccDist=false; plot_range_vs_control
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; diary off; more on; pack
persub = @ts_jas; doWaves=false; ctlvar={'ngdc_depth',@nanmean,9,9,10}; doAccDist=false; plot_range_vs_control
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; diary off; more on; pack
persub = @ts_jfm; doWaves=false; ctlvar={'ngdc_depth',@nanmean,9,9,10}; doAccDist=true; plot_range_vs_control
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; diary off; more on; pack
persub = @ts_jas; doWaves=false; ctlvar={'ngdc_depth',@nanmean,9,9,10}; doAccDist=true; plot_range_vs_control
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; diary off; more on; pack

persub = @ts_jfm; HsRange=[0,0.5]; ctlvar={'ngdc_depth',@nanmean,9,9,10}; doAccDist=false; plot_range_vs_control
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; diary off; more on; pack
persub = @ts_jas; HsRange=[0,0.5]; ctlvar={'ngdc_depth',@nanmean,9,9,10}; doAccDist=false; plot_range_vs_control
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; diary off; more on; pack
persub = @ts_jfm; HsRange=[0,0.5]; ctlvar={'ngdc_depth',@nanmean,9,9,10}; doAccDist=true; plot_range_vs_control
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; diary off; more on; pack
persub = @ts_jas; HsRange=[0,0.5]; ctlvar={'ngdc_depth',@nanmean,9,9,10}; doAccDist=true; plot_range_vs_control
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; diary off; more on; pack


%% SEAFLOOR SLOPE
persub = @ts_jfm; doWaves=false; ctlvar={'ngdc_beta',@nanmean,9,9,10}; doAccDist=false; plot_range_vs_control
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; diary off; more on; pack
persub = @ts_jas; doWaves=false; ctlvar={'ngdc_beta',@nanmean,9,9,10}; doAccDist=false; plot_range_vs_control
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; diary off; more on; pack
persub = @ts_jfm; doWaves=false; ctlvar={'ngdc_beta',@nanmean,9,9,10}; doAccDist=true; plot_range_vs_control
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; diary off; more on; pack
persub = @ts_jas; doWaves=false; ctlvar={'ngdc_beta',@nanmean,9,9,10}; doAccDist=true; plot_range_vs_control
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; diary off; more on; pack

persub = @ts_jfm; HsRange=[0,0.5]; ctlvar={'ngdc_beta',@nanmean,9,9,10}; doAccDist=false; plot_range_vs_control
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; diary off; more on; pack
persub = @ts_jas; HsRange=[0,0.5]; ctlvar={'ngdc_beta',@nanmean,9,9,10}; doAccDist=false; plot_range_vs_control
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; diary off; more on; pack
persub = @ts_jfm; HsRange=[0,0.5]; ctlvar={'ngdc_beta',@nanmean,9,9,10}; doAccDist=true; plot_range_vs_control
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; diary off; more on; pack
persub = @ts_jas; HsRange=[0,0.5]; ctlvar={'ngdc_beta',@nanmean,9,9,10}; doAccDist=true; plot_range_vs_control
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; diary off; more on; pack
