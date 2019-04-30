1;
% MY_STARTUP script
% 
% Normally called from MATLABROOT/toolbox/local/startup.m
% (Or formerly from MATLABROOT/work/startup.m)
% 
% dedalus@alum.mit.edu, lgramer@upstreamtechnology.com, lew.gramer@noaa.gov, 2001-2017
% 
% Last Saved Time-stamp: <Sat 2017-12-09 14:32:51 Eastern Standard Time gramer>


%
% Options
%

format compact;
format short;
%more on;
more off;


%
% Path modifications
%

%% TITAN XP
%MYDOCS = 'C:\Documents and Settings\lew.gramer\My Documents\';
%% TETHYS Win7
%MYDOCS = 'C:\Users\gramer\Documents\';
%% TITAN Win7
%MYDOCS = 'C:\Users\lew.gramer\Documents\';
%MYDATA = 'E:\';

% MANANNAN Win7
MYDOCS = 'C:\Users\gramer\Documents\';
MYDATA = 'D:\';

%% MUIRGEN macOS
%MYDOCS = '~/';
%MYDATA = '~/';

MHOME = [MYDOCS 'MATLAB\'];

%%path(path,[MYDOCS 'RSMAS\Coastal\thesis\src\']);
%path(path,[MYDATA 'thesis\src\']);
path(path,[MYDATA 'heat_budget\src\']);

%path(path,[MYDOCS 'coral\CRCP\Sediment\']);
path(path,[MYDATA 'coral\']);
noaa_aoml_coral_project_path;

% All code should migrate to Ecoforecasts or Heat_Budget toolkits
%%%% path(path,[MYDATA 'Postdoc\']);

% Are we running in Octave (free software) vs. MATLAB (*not*)?
if ( ~exist('isOctave') )
  disp('Assuming this is NOT Octave');
  isOctave = false;
end;
if ( isOctave )
  % pkg load auto;
  % pkg load data-smoothing
  % pkg load dataframe
  % pkg load financial
  % pkg load fl-core
  % pkg load fuzzy-logic-toolkit
  % pkg load image
  % pkg load io
  % pkg load linear-algebra
  % pkg load lssa % Spectral decompositions of irregularly-spaced time series
  % pkg load ltfat % Large Time/Frequency Analysis Toolbox
  pkg load netcdf
  pkg load signal
  % pkg load splines
  pkg load statistics
  % pkg load stk % (not so) Small Toolbox for Kriging
  % pkg load strings
  % pkg load struct
  % pkg load zenity % Simple graphical user interfaces
else
  if ( exist('opengl') > 1 )
    %opengl hardware
    %opengl hardwarebasic
    opengl software
  end;
end;

% Carnegie-Mellon University tools (units, color names in plotting, finite difference derivative functions, and easy installation of Matlab functions)
%path(path,[MHOME '+cmu/']);
import cmu.*

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
path(path,[MYDATA,'ecoforecasts/']);
% Seafloor topography
path(path,[MYDATA,'ecoforecasts/coast/']);
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
if ( ~isOctave )
  % JLAB: Matlab freeware for data analysis: http://www.jmlilly.net/jmlsoft.html
  path(path,[MHOME,'jlab/']);
end;
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

% The CSIRO netCDF/OPeNDAP interface to MATLAB:
%   http://www.marine.csiro.au/sw/matlab-netcdf.html
% CSIROTBX='matlab_netCDF_OPeNDAP\';
% javaaddpath([MHOME CSIROTBX 'toolsUI-4.0.49.jar'],'-end');
% path(path,[MHOME,CSIROTBX]);
CSIROTBX='matlab_netCDF_OPeNDAP\';
if ( isOctave )
  javaaddpath([MHOME CSIROTBX 'toolsUI-4.0.49.jar']);
  %javaaddpath([MHOME CSIROTBX 'toolsUI-4.6.jar']);
else
  javaaddpath([MHOME CSIROTBX 'toolsUI-4.0.49.jar'],'-end');
  %javaaddpath([MHOME CSIROTBX 'toolsUI-4.6.jar'],'-end');
end;
path(path,[MHOME,CSIROTBX]);
clear CSIROTBX

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
if ( isOctave )
  javaaddpath([MHOME NJTBX 'toolsUI-4.0.49.jar']);
  javaaddpath([MHOME NJTBX 'njTools-2.0.12_jre1.5.jar']);
  % javaaddpath([MHOME NJTBX 'njTools-2.0.12_jre1.6.jar']);
else
  javaaddpath([MHOME NJTBX 'toolsUI-4.0.49.jar'],'-end');
  javaaddpath([MHOME NJTBX 'njTools-2.0.12_jre1.5.jar'],'-end');
  % javaaddpath([MHOME NJTBX 'njTools-2.0.12_jre1.6.jar'],'-end');

  % logger = org.apache.log4j.Logger.getLogger('org.apache.fop');
  % logger.addAppender(org.apache.log4j.varia.NullAppender);
  % logger = org.apache.log4j.Logger.getLogger('org.apache.xmlgraphics.image');
  % logger.addAppender(org.apache.log4j.varia.NullAppender);
  % logger = org.apache.log4j.Logger.getLogger('FOP');
  % logger.addAppender(org.apache.log4j.varia.NullAppender);
end;
path(path,[MHOME,NJTBX]);
path(path,[MHOME,NJTBX,'examples/']);
path(path,[MHOME,NJTBX,'njFunc/']);
path(path,[MHOME,NJTBX,'njTBX-2.0/']);
path(path,[MHOME,NJTBX,'njTBX-2.0/Utilities/']);
clear NJTBX

% Nortek data processing tools for Aquadopp, AWAC, etc.
path(path,[MHOME,'Nortek/']);
% Pattern Recognition Tools from http://www.prtools.org
warning('OFF','MATLAB:dispatcher:nameConflict');
path(path,[MHOME,'prtools/']);
warning('ON','MATLAB:dispatcher:nameConflict');
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
%initdirs(MHOME);
initdirs;

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

clear MHOME
