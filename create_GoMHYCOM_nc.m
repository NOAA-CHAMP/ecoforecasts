function flds = read_hycom_gom_reanalysis()
%function flds = read_hycom_gom_reanalysis()
% CREATE dl-archive/GoMHYCOM.nc file with latitude and longitude from DODS (expt_32.5)

  datapath = get_ecoforecasts_path('dl-archive');

  % RESULTS FROM: http://tds.hycom.org/thredds/dodsC/GOMl0.04/expt_32.5/hrly.ascii?Latitude[0:1:384],Longitude[0:1:540]
  lats = [ ...
      18.091648, 18.129667, 18.167677, 18.205679, 18.243671, 18.281658, 18.319633, 18.357603, 18.395563, 18.433516, 18.471458, 18.509394, 18.54732, 18.585238, 18.623148, 18.661049, 18.698942, 18.736826, 18.774702, 18.81257, 18.85043, 18.888279, 18.92612, 18.963955, 19.00178, 19.039595, 19.077402, 19.115202, 19.15299, 19.190773, 19.228546, 19.26631, 19.304066, 19.341812, 19.379549, 19.417278, 19.455, 19.492712, 19.530415, 19.56811, 19.605795, 19.64347, 19.681139, 19.718798, 19.756447, 19.794088, 19.83172, 19.869343, 19.906958, 19.944563, 19.98216, 20.019747, 20.057325, 20.094894, 20.132456, 20.170008, 20.20755, 20.245083, 20.282608, 20.320122, 20.357628, 20.395126, 20.432613, 20.470093, 20.507563, 20.545023, 20.582474, 20.619915, 20.657349, 20.694773, 20.732187, 20.769592, 20.806988, 20.844374, 20.881752, 20.91912, 20.956478, 20.993828, 21.031168, 21.068499, 21.10582, 21.143133, 21.180435, 21.217728, 21.255013, 21.292286, 21.329552, 21.366806, 21.404053, 21.441288, 21.478516, 21.515734, 21.55294, 21.59014, 21.627329, 21.664507, 21.701677, 21.738836, 21.775988, 21.81313, 21.85026, 21.88738, 21.924494, 21.961596, 21.998688, 22.03577, 22.072844, 22.109907, 22.146961, 22.184006, 22.221039, 22.258064, 22.295078, 22.332083, 22.369078, 22.406063, 22.443039, 22.480003, 22.51696, 22.553905, 22.590841, 22.627768, 22.664682, 22.70159, 22.738485, 22.775372, 22.812248, 22.849113, 22.88597, 22.922815, 22.959652, 22.996479, 23.033295, 23.0701, 23.106897, 23.143682, 23.180458, 23.217224, 23.25398, 23.290726, 23.327461, 23.364185, 23.4009, 23.437605, 23.4743, 23.510984, 23.547659, 23.584324, 23.620977, 23.65762, 23.694254, 23.730877, 23.767488, 23.804092, 23.840683, 23.877266, 23.913837, 23.950397, 23.986948, 24.023489, 24.060019, 24.096539, 24.133047, 24.169546, 24.206034, 24.242512, 24.27898, 24.315437, 24.351883, 24.388319, 24.424746, 24.461159, 24.497564, 24.533958, 24.570341, 24.606714, 24.643076, 24.679428, 24.715769, 24.7521, 24.78842, 24.824728, 24.861027, 24.897314, 24.933592, 24.969858, 25.006115, 25.04236, 25.078594, 25.114817, 25.151031, 25.187233, 25.223425, 25.259605, 25.295774, 25.331934, 25.368082, 25.40422, 25.440348, 25.476463, 25.512568, 25.548662, 25.584745, 25.620817, 25.65688, 25.69293, 25.72897, 25.765, 25.801016, 25.837023, 25.87302, 25.909004, 25.944979, 25.980942, 26.016893, 26.052835, 26.088764, 26.124683, 26.160591, 26.19649, 26.232374, 26.26825, 26.304113, 26.339966, 26.375807, 26.411638, 26.447458, 26.483265, 26.519062, 26.554848, 26.590624, 26.626387, 26.662138, 26.69788, 26.73361, 26.76933, 26.805037, 26.840733, 26.876417, 26.91209, 26.947754, 26.983404, 27.019045, 27.054672, 27.09029, 27.125896, 27.161491, 27.197075, 27.232647, 27.268206, 27.303755, 27.339294, 27.37482, 27.410336, 27.445839, 27.48133, 27.516811, 27.55228, 27.587738, 27.623184, 27.658619, 27.694044, 27.729456, 27.764856, 27.800245, 27.835623, 27.870987, 27.906342, 27.941685, 27.977016, 28.012335, 28.047644, 28.082941, 28.118225, 28.153498, 28.18876, 28.22401, 28.259249, 28.294476, 28.32969, 28.364893, 28.400085, 28.435266, 28.470434, 28.50559, 28.540735, 28.575869, 28.61099, 28.646101, 28.681198, 28.716284, 28.751358, 28.78642, 28.821472, 28.856512, 28.891539, 28.926554, 28.96156, 28.99655, 29.03153, 29.066498, 29.101456, 29.136398, 29.171331, 29.206253, 29.241161, 29.276058, 29.310944, 29.345816, 29.380678, 29.415527, 29.450363, 29.48519, 29.520002, 29.554804, 29.589594, 29.62437, 29.659136, 29.69389, 29.72863, 29.763361, 29.798077, 29.832783, 29.867476, 29.902157, 29.936827, 29.971483, 30.006128, 30.04076, 30.075382, 30.109991, 30.144587, 30.17917, 30.213743, 30.248304, 30.282852, 30.317387, 30.351912, 30.386423, 30.420921, 30.455408, 30.489883, 30.524345, 30.558796, 30.593235, 30.62766, 30.662075, 30.696476, 30.730865, 30.765242, 30.799606, 30.83396, 30.8683, 30.902626, 30.936943, 30.971245, 31.005537, 31.039816, 31.074081, 31.108335, 31.142576, 31.176805, 31.211023, 31.245228, 31.279419, 31.313599, 31.347765, 31.38192, 31.416063, 31.450193, 31.48431, 31.518415, 31.552507, 31.586588, 31.620657, 31.65471, 31.688755, 31.722784, 31.756802, 31.790808, 31.8248, 31.858782, 31.892748, 31.926704, 31.960648 ...
         ];
  lons = [ ...
      -98.0, -97.96002, -97.91998, -97.880005, -97.839966, -97.79999, -97.76001, -97.71997, -97.67999, -97.640015, -97.599976, -97.56, -97.52002, -97.47998, -97.44, -97.400024, -97.359985, -97.32001, -97.28003, -97.23999, -97.20001, -97.160034, -97.119995, -97.08002, -97.03998, -97.0, -96.96002, -96.91998, -96.880005, -96.839966, -96.79999, -96.76001, -96.71997, -96.67999, -96.640015, -96.599976, -96.56, -96.52002, -96.47998, -96.44, -96.400024, -96.359985, -96.32001, -96.28003, -96.23999, -96.20001, -96.160034, -96.119995, -96.08002, -96.03998, -96.0, -95.96002, -95.91998, -95.880005, -95.839966, -95.79999, -95.76001, -95.71997, -95.67999, -95.640015, -95.599976, -95.56, -95.52002, -95.47998, -95.44, -95.400024, -95.359985, -95.32001, -95.28003, -95.23999, -95.20001, -95.160034, -95.119995, -95.08002, -95.03998, -95.0, -94.96002, -94.91998, -94.880005, -94.839966, -94.79999, -94.76001, -94.71997, -94.67999, -94.640015, -94.599976, -94.56, -94.52002, -94.47998, -94.44, -94.400024, -94.359985, -94.32001, -94.28003, -94.23999, -94.20001, -94.160034, -94.119995, -94.08002, -94.03998, -94.0, -93.96002, -93.91998, -93.880005, -93.839966, -93.79999, -93.76001, -93.71997, -93.67999, -93.640015, -93.599976, -93.56, -93.52002, -93.47998, -93.44, -93.400024, -93.359985, -93.32001, -93.28003, -93.23999, -93.20001, -93.160034, -93.119995, -93.08002, -93.03998, -93.0, -92.96002, -92.91998, -92.880005, -92.839966, -92.79999, -92.76001, -92.71997, -92.67999, -92.640015, -92.599976, -92.56, -92.52002, -92.47998, -92.44, -92.400024, -92.359985, -92.32001, -92.28003, -92.23999, -92.20001, -92.160034, -92.119995, -92.08002, -92.03998, -92.0, -91.96002, -91.91998, -91.880005, -91.839966, -91.79999, -91.76001, -91.71997, -91.67999, -91.640015, -91.599976, -91.56, -91.52002, -91.47998, -91.44, -91.400024, -91.359985, -91.32001, -91.28003, -91.23999, -91.20001, -91.160034, -91.119995, -91.08002, -91.03998, -91.0, -90.96002, -90.91998, -90.880005, -90.839966, -90.79999, -90.76001, -90.71997, -90.67999, -90.640015, -90.599976, -90.56, -90.52002, -90.47998, -90.44, -90.400024, -90.359985, -90.32001, -90.28003, -90.23999, -90.20001, -90.160034, -90.119995, -90.08002, -90.03998, -90.0, -89.96002, -89.91998, -89.880005, -89.839966, -89.79999, -89.76001, -89.71997, -89.67999, -89.640015, -89.599976, -89.56, -89.52002, -89.47998, -89.44, -89.400024, -89.359985, -89.32001, -89.28003, -89.23999, -89.20001, -89.160034, -89.119995, -89.08002, -89.03998, -89.0, -88.96002, -88.91998, -88.880005, -88.839966, -88.79999, -88.76001, -88.71997, -88.67999, -88.640015, -88.599976, -88.56, -88.52002, -88.47998, -88.44, -88.400024, -88.359985, -88.32001, -88.28003, -88.23999, -88.20001, -88.160034, -88.119995, -88.08002, -88.03998, -88.0, -87.96002, -87.91998, -87.880005, -87.839966, -87.79999, -87.76001, -87.71997, -87.67999, -87.640015, -87.599976, -87.56, -87.52002, -87.47998, -87.44, -87.400024, -87.359985, -87.32001, -87.28003, -87.23999, -87.20001, -87.160034, -87.119995, -87.08002, -87.03998, -87.0, -86.96002, -86.91998, -86.880005, -86.839966, -86.79999, -86.76001, -86.71997, -86.67999, -86.640015, -86.599976, -86.56, -86.52002, -86.47998, -86.44, -86.400024, -86.359985, -86.32001, -86.28003, -86.23999, -86.20001, -86.160034, -86.119995, -86.08002, -86.03998, -86.0, -85.96002, -85.91998, -85.880005, -85.839966, -85.79999, -85.76001, -85.71997, -85.67999, -85.640015, -85.599976, -85.56, -85.52002, -85.47998, -85.44, -85.400024, -85.359985, -85.32001, -85.28003, -85.23999, -85.20001, -85.160034, -85.119995, -85.08002, -85.03998, -85.0, -84.96002, -84.91998, -84.880005, -84.839966, -84.79999, -84.76001, -84.71997, -84.67999, -84.640015, -84.599976, -84.56, -84.52002, -84.47998, -84.44, -84.400024, -84.359985, -84.32001, -84.28003, -84.23999, -84.20001, -84.160034, -84.119995, -84.08002, -84.03998, -84.0, -83.96002, -83.91998, -83.880005, -83.839966, -83.79999, -83.76001, -83.71997, -83.67999, -83.640015, -83.599976, -83.56, -83.52002, -83.47998, -83.44, -83.400024, -83.359985, -83.32001, -83.28003, -83.23999, -83.20001, -83.160034, -83.119995, -83.08002, -83.03998, -83.0, -82.96002, -82.91998, -82.880005, -82.839966, -82.79999, -82.76001, -82.71997, -82.67999, -82.640015, -82.599976, -82.56, -82.52002, -82.47998, -82.44, -82.400024, -82.359985, -82.32001, -82.28003, -82.23999, -82.20001, -82.160034, -82.119995, -82.08002, -82.03998, -82.0, -81.96002, -81.91998, -81.880005, -81.839966, -81.79999, -81.76001, -81.71997, -81.67999, -81.640015, -81.599976, -81.56, -81.52002, -81.47998, -81.44, -81.400024, -81.359985, -81.32001, -81.28003, -81.23999, -81.20001, -81.160034, -81.119995, -81.08002, -81.03998, -81.0, -80.96002, -80.91998, -80.880005, -80.839966, -80.79999, -80.76001, -80.71997, -80.67999, -80.640015, -80.599976, -80.56, -80.52002, -80.47998, -80.44, -80.400024, -80.359985, -80.32001, -80.28003, -80.23999, -80.20001, -80.160034, -80.119995, -80.08002, -80.03998, -80.0, -79.96002, -79.91998, -79.880005, -79.839966, -79.79999, -79.76001, -79.71997, -79.67999, -79.640015, -79.599976, -79.56, -79.52002, -79.47998, -79.44, -79.400024, -79.359985, -79.32001, -79.28003, -79.23999, -79.20001, -79.160034, -79.119995, -79.08002, -79.03998, -79.0, -78.96002, -78.91998, -78.880005, -78.839966, -78.79999, -78.76001, -78.71997, -78.67999, -78.640015, -78.599976, -78.56, -78.52002, -78.47998, -78.44, -78.400024, -78.359985, -78.32001, -78.28003, -78.23999, -78.20001, -78.160034, -78.119995, -78.08002, -78.03998, -78.0, -77.96002, -77.91998, -77.880005, -77.839966, -77.79999, -77.76001, -77.71997, -77.67999, -77.640015, -77.599976, -77.56, -77.52002, -77.47998, -77.44, -77.400024, -77.359985, -77.32001, -77.28003, -77.23999, -77.20001, -77.160034, -77.119995, -77.08002, -77.03998, -77.0, -76.96002, -76.91998, -76.880005, -76.839966, -76.79999, -76.76001, -76.71997, -76.67999, -76.640015, -76.599976, -76.56, -76.52002, -76.47998, -76.44, -76.400024 ...
         ];

  fname = fullfile(datapath,'GoMHYCOM.nc');

  %netcdf.open(fname,'WRITE');
  %netcdf.reDef(ncid);
  ncid = netcdf.create(fname,'CLOBBER');

  latDimID = netcdf.defDim(ncid,'latitude',numel(lats));
  lonDimID = netcdf.defDim(ncid,'longitude',numel(lons));

  latVarID = netcdf.defVar(ncid,'latitude','double',latDimID);
  lonVarID = netcdf.defVar(ncid,'longitude','double',lonDimID);

  netcdf.endDef(ncid);
  netcdf.putVar(ncid,latVarID,lats);
  netcdf.putVar(ncid,lonVarID,lons);

  netcdf.close(ncid);

return;
