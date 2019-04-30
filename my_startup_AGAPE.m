% 
% Normally called from MATLABROOT/work/startup.m
% 
% 
% Last Saved Time-stamp: <Fri 2015-02-13 16:43:09 Eastern Standard Time gramer>


%
% Options
%

format compact;
format short;
more on;


%
% Path modifications
%

%% TITAN XP
%MYDOCS = 'C:\Documents and Settings\lew.gramer\My Documents\';
% TITAN Win7
MYDOCS = 'C:\Users\gramer\Documents\';

MHOME = [MYDOCS 'MATLAB\'];

path(path,[MYDOCS 'Postdoc\']);

path(path,[MYDOCS 'RSMAS\Coastal\thesis\src\']);


% Air-Sea Interaction Toolbox
path(path,[MHOME 'air_sea/']);
% Fourier analysis tools developed at RSMAS AMP - sent me by Xiaofang Zhu
path(path,[MHOME 'AMPtools/']);
% Bootstrap Toolbox by AM Zoubir & DR Iskander, Curtin University of Technology-Australia
path(path,[MHOME,'bootstrap/']);
% % Lew Gramer's tools for creating, analyzing, and plotting HydroBase2
% % (isopycnal-averaged) and Levitus (1-deg, 1/4-deg) ocean climatologies
% path(path,[MHOME,'climatology/']);
% % Derrick Snowden's unique COLORMAPs for oceanographic variables
% path(path,[MHOME,'colormaps/']);
% Like UNIQUE(), but more powerful and with fuzzy tolerance
path(path,[MHOME,'consolidator/']);
% Lew's ICON/MATLAB Ecoforecasting Toolkit
path(path,[MHOME 'ecoforecasts/']);
% Seafloor topography
path(path,[MHOME 'ecoforecasts/coast/']);
% % LeSage Econometrics Toolbox: GARCH, Bayesian Markov Chain Monte Carlo,
% % Gibbs sampling, simplex optimization, et al: http://www.econ.utoledo.edu
% path(path,[MHOME 'Econometrics/']);
% Regular, rotated, reduced, teleconnection decomposed, complex, and "other" EOFs
% by Martijn Hooimeijer (1998) http://hydr.ct.tudelft.nl/wbk/public/hooimeijer/
path(path,[MHOME 'EOFs/']);
% TOGA-COARE Heat Flux calculations
path(path,[MHOME 'fairall/']);
% Fairall Radiative Flux calculations: ftp://ftp1.esrl.noaa.gov/users/cfairall/radflux/matlab/
path(path,[MHOME 'fairall/radflux/']);
% Fast Independent Component Analysis
path(path,[MHOME 'FastICA/']);
% Gibbs Sea Water Oceanographic Toolbox: http://www.TEOS-10.org
% McDougall T. J. and P. M. Barker, 2011: Getting started with TEOS-10 and the Gibbs Seawater (GSW) Oceanographic Toolbox, 28pp., SCOR/IAPSO WG127, ISBN 978-0-646-55621-5.
path(path,[MHOME,'GSW/']);
path(path,[MHOME,'GSW/html']);
path(path,[MHOME,'GSW/library']);
path(path,[MHOME,'GSW/pdf']);
% Growing Hierarchical SOM (see somboolbox below)
path(path,[MHOME,'ghsom/']);
% Hilbert-Huang Transform - signal decomposition and marginal spectra
path(path,[MHOME,'hilberthuang/']);
% Independent Component Analysis, Sammon nets, other stuff
path(path,[MHOME,'icasso122/']);
% JLAB: Matlab freeware for data analysis: http://www.jmlilly.net/jmlsoft.html
path(path,[MHOME,'jlab/']);
% CSIRO's IMOS-Toolbox from: http://code.google.com/p/imos-toolbox/source/browse/#svn%2Ftrunk
path(path,[MHOME,'imos-toolbox/']);
%path(path,[MHOME,'imos-toolbox/AutomaticQC/']);
%path(path,[MHOME,'imos-toolbox/DDB/']);
%path(path,[MHOME,'imos-toolbox/FlowManager/']);
%path(path,[MHOME,'imos-toolbox/GUI/']);
%path(path,[MHOME,'imos-toolbox/Graph/']);
% For 3.2 and after?
path(path,[MHOME,'imos-toolbox/IMOS/']);
%path(path,[MHOME,'imos-toolbox/Java/']);
%path(path,[MHOME,'imos-toolbox/NetCDF/']);
path(path,[MHOME,'imos-toolbox/Parser/']);
%path(path,[MHOME,'imos-toolbox/Preprocessing/']);
%path(path,[MHOME,'imos-toolbox/Seawater/']);
%path(path,[MHOME,'imos-toolbox/Tests/']);
path(path,[MHOME,'imos-toolbox/Util/']);
% % Matlab Structs Tool (LOADDAP): http://opendap.org/download/ml-structs.html#3.6.1
% path(path,'c:/opendap/loaddap/');
% % Adaptive low-pass filtering, v. Mann (2008) Smoothing of climate time series revisited.
% path(path,[MHOME 'Mann_filtering/']);
% Various image processing/Computer Vision tools
path(path,[MHOME 'MatlabFns/Match/']);
path(path,[MHOME,'netlab/']);
%
% % NetCDF Java Toolbox (OLD)
% javaaddpath([MHOME 'njToolbox-2.0\toolsUI-4.0.49.jar'],'-end');
% javaaddpath([MHOME 'njToolbox-2.0\njTools-2.0.12_jre1.5.jar'],'-end');
% % javaaddpath([MHOME 'njToolbox-2.0\njTools-2.0.12_jre1.6.jar'],'-end');
% path(path,[MHOME,'njToolbox-2.0/']);
% path(path,[MHOME,'njToolbox-2.0/examples/']);
% path(path,[MHOME,'njToolbox-2.0/njFunc/']);
% path(path,[MHOME,'njToolbox-2.0/njTBX-2.0/']);
% path(path,[MHOME,'njToolbox-2.0/njTBX-2.0/Utilities/']);
% %
% NetCDF Java Toolbox: http://sourceforge.net/apps/trac/njtbx
NJTBX='matlab-njTbx-2.0.05/';
javaaddpath([MHOME NJTBX 'toolsUI-4.0.49.jar'],'-end');
javaaddpath([MHOME NJTBX 'njTools-2.0.12_jre1.5.jar'],'-end');
% javaaddpath([MHOME NJTBX 'njTools-2.0.12_jre1.6.jar'],'-end');
path(path,[MHOME,NJTBX]);
path(path,[MHOME,NJTBX,'examples/']);
path(path,[MHOME,NJTBX,'njFunc/']);
path(path,[MHOME,NJTBX,'njTBX-2.0/']);
path(path,[MHOME,NJTBX,'njTBX-2.0/Utilities/']);
clear NJTBX
%
% Nortek data processing tools for Aquadopp, AWAC, etc.
path(path,[MHOME,'Nortek/']);
% Pattern Recognition Tools from http://www.prtools.org
path(path,[MHOME,'prtools/']);
path(path,[MHOME,'oceans/']);
% WMO GRiB-format file reader (via NOAA NCEP WaveWatch III site:
%   http://polar.ncep.noaa.gov/waves/download.shtml
path(path,[MHOME,'read_grib/']);
% Smoothened Data-Histogram (used on SOM, see below)
path(path,[MHOME,'sdh/']);
path(path,[MHOME,'snackbar/']);
% % Wiberg & Sherwood 2008: Calculating wave-generated bottom orbital velocities from surface-wave parameters
% path(path,[MHOME,'snackbar/wiberg_wave_bottom_vel/']);
% Self-Organizing (or Kohonen) Maps toolbox
path(path,[MHOME,'somtoolbox/']);
% Sea-Water toolbox
path(path,[MHOME,'Sw/']);
% RPS T_TIDE: http://www.eos.ubc.ca/~rich/
% R. Pawlowicz, B. Beardsley, and S. Lentz, "Classical tidal harmonic analysis
% including error estimates in MATLAB using T_TIDE", Computers and Geosciences
% 28 (2002), 929-937.
path(path,[MHOME,'t_tide/']);
% Tide Model Driver toolbox from: http://volkov.oce.orst.edu/tides/
path(path,[MHOME,'tmd_toolbox/']);
% UTide - Unified Tidal Analysis and Prediction: http://www.po.gso.uri.edu/~codiga/utide/utide.htm
path(path,[MHOME,'UTide/']);
% Grinsted implementation - Wavelet analysis
path(path,[MHOME 'wavelets/']);
path(path,[MHOME 'wavelets/wtc/']);
path(path,[MHOME 'wavelets/wavelet/']);

clear MHOME



% PHOME = '\\cygnus\gramer\home\matlab\';

% path(path,[PHOME]);

% % Like UNIQUE(), but more powerful and with fuzzy tolerance
% path(path,[PHOME,'consolidator/']);
% path(path,[PHOME,'ecoforecasts/']);
% path(path,[PHOME,'ecoforecasts/coast/']);
% path(path,[PHOME 'EOFs/']);


% path(path,[PHOME,'ADCP/']);
% path(path,[PHOME,'Advstats']);
% path(path,[PHOME,'AirSea/']);
% path(path,[PHOME,'Arfit/']);
% path(path,[PHOME,'Bsplines']);
% path(path,[PHOME,'Claudia/']);
% path(path,[PHOME,'climatology']);
% path(path,[PHOME,'CODAS3/']);
% path(path,[PHOME,'compstats']);
% path(path,[PHOME,'Colormap']);
% path(path,[PHOME,'coral/']);
% path(path,[PHOME,'CTDCalib/']);
% path(path,[PHOME,'DataSet/']);
% path(path,[PHOME,'DateTime']);
% path(path,[PHOME,'Dsatbx/']);
% path(path,[PHOME,'Dvt/']);
% path(path,[PHOME,'DynModes/']);
% path(path,[PHOME,'EarthSci/']);
% path(path,[PHOME,'Epic/']);
% path(path,[PHOME,'Eps/']);
% path(path,[PHOME,'Even/']);
% path(path,[PHOME,'Fmex/']);
% path(path,[PHOME,'Geodetics/']);
% path(path,[PHOME,'Geography']);
% path(path,[PHOME,'gkslib/']);
% path(path,[PHOME,'Glmlab/']);
% path(path,[PHOME,'Graphics']);
% path(path,[PHOME,'Hist']);
% path(path,[PHOME,'HtmlTool']);
% path(path,[PHOME,'ifmObana/']);
% path(path,[PHOME,'Imputation/']);
% path(path,[PHOME,'Integration']);
% path(path,[PHOME,'IO/']);
% path(path,[PHOME,'MClasses']);
% path(path,[PHOME,'MClasses/Generic']);
% path(path,[PHOME,'MClasses/Interval']);
% path(path,[PHOME,'MClasses/Matfile']);
% path(path,[PHOME,'MClasses/Mode']);
% path(path,[PHOME,'MClasses/Mooring']);
% path(path,[PHOME,'MClasses/Vos']);
% path(path,[PHOME,'MClasses/batch']);
% path(path,[PHOME,'MClasses/presto']);
% path(path,[PHOME,'MClasses/time_series']);
% path(path,[PHOME,'MClasses/Vivace/']);
% path(path,[PHOME,'MClasses/Profile']);
% path(path,[PHOME,'Makehelp/']);
% path(path,[PHOME,'mat2html/']);
% path(path,[PHOME,'matlab_netcdf_5_0/']);
% path(path,[PHOME,'mexcdf60/']);
% path(path,[PHOME,'MexEPS/']);

% %%%% ??? MAYBE NEEDED AFTER ALL
% path(path,[PHOME,'mexnc']);

% path(path,[PHOME,'Misc/']);
% path(path,[PHOME,'m_map/']);
% path(path,[PHOME,'m_map/m_namebox/']);
% path(path,[PHOME,'MySQL']);
% path(path,[PHOME,'netcdf']);
% path(path,[PHOME,'netcdf/ncfiles']);
% path(path,[PHOME,'netcdf/nctype']);
% path(path,[PHOME,'netcdf/ncutility']);

% %%%% ??? MAYBE NEEDED AFTER ALL
% path(path,[PHOME,'netcdf_toolbox']);
% path(path,[PHOME,'netcdf_toolbox/netcdf']);
% % NOTE: This path differs from that mentioned in README... (LGramer, 2005-06-29)
% path(path,[PHOME,'netcdf_toolbox/netcdf/ncsource']);
% path(path,[PHOME,'netcdf_toolbox/netcdf/nctype']);
% path(path,[PHOME,'netcdf_toolbox/netcdf/ncutility']);

% path(path,[PHOME,'NumMethods/']);
% path(path,[PHOME,'Nurbs/']);
% path(path,[PHOME,'Oceanography']);
% path(path,[PHOME,'Oceanography/tsplot']);
% path(path,[PHOME,'Oceanography/Hydrobase']);
% path(path,[PHOME,'OMP2/']);
% path(path,[PHOME,'omviz/']);
% path(path,[PHOME,'OPNML/']);
% path(path,[PHOME,'PCCruiseCD/']);
% path(path,[PHOME,'RDADCP/']);
% path(path,[PHOME,'Regularization/']);
% path(path,[PHOME,'R12/']);
% path(path,[PHOME,'Smooth/']);
% path(path,[PHOME,'snackbar']);
% % Old tools for netCDF reading and creation
% path(path,[PHOME,'snctools']);
% path([PHOME,'Statistics'],path);
% path(path,[PHOME,'Strings']);
% path(path,[PHOME,'SaGA/']);
% path(path,[PHOME,'Statbx40']);
% path(path,[PHOME,'Spatial']);
% path(path,[PHOME,'Splines/']);
% path(path,[PHOME,'Staplot/']);
% path(path,[PHOME,'Steger/']);

%%%% ??? NOTE: This 'Sw' is older, but includes additional components from
%%%% ??? NOTE:  URI for sound-velocity and deep-ocean density calculations.
% path(path,[PHOME,'Sw/']);

% path(path,[PHOME,'timeplt5/']);
% path(path,[PHOME,'TimeSeries/']);
% path(path,[PHOME,'Todd/']);
% path(path,[PHOME,'transports']);
% path(path,[PHOME,'WaveCov/']);
% path(path,[PHOME,'WaveLab802/']);
% path(path,[PHOME,'Wavelets/']);
% path(path,[PHOME,'wetcdf/']);
% path(path,[PHOME,'XmlStuff/']);

% addutils;
%%%% ??? NOTE: Commenting out ADDUTILS above actually keeps the following
%%%% ??? NOTE: *TWELVE* utility directories from being added to the path.
% datautil
% fileutil
% graphutil
% imgutil
% mathutil
% matutil
% numutil
% polyutil
% statutil
% strutil
% sysutil
% timeutil

clear PHOME


% %%%% ??? MAYBE NEEDED AFTER ALL
% %
% % netCDF initialization
% %
% global nctbx_options;
% nctbx_options.theAutoscale = 1;
% nctbx_options.theAutoNaN = 1;


% %
% % Keep tabs on the current path on this host...
% %
% fid = fopen('~/.matlab.path', 'w+');
% if ( fid > 0 )
%     fprintf(fid, '%s\n', path);
%     fclose(fid);
% end
% clear fid;


% Initialize the directory stack
initdirs;
