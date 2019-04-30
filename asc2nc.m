function asc2nc(ascfname)

  error('Function not ready for prime time yet');

  ncfname = strrep(ascfname,'.asc','.nc');
  nccreate(ncfname, 'bathymetry', 'Dimensions', {'lon' 18 'lat' 18}, 'Format','classic')

return;
