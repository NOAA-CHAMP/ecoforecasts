1;
%% SCRIPT read_coral_and_hard_bottom.m:
% Read SHAPE files with polygons outlining all coral and hard-bottom habitat
% from the FWC Florida Unified Coral Reef Habitat Map v2.0


coastpath = get_ecoforecasts_path('coast');

matfname = fullfile(coastpath,['coral_and_hard_bottom_polygons.mat']);

if ( exist(matfname,'file') )
  disp(['Loading ',matfname]);
  load(matfname);

else
  disp('Reading SHAPE file');
  shps = shaperead(fullfile(coastpath,'Coral_and_Hard_Bottom_Habitats_in_Florida','Coral_and_Hard_Bottom_Habitats_in_Florida'));

  rfhd.lon = []; rfhd.lat = [];
  reef.lon = []; reef.lat = [];
  rfhd.polylon = []; rfhd.polylat = [];
  reef.polylon = []; reef.polylat = [];

  disp('Parsing SHAPE structure');
  for six=1:numel(shps)
    n = numel(shps(six).X);
    rfhd.lon(end+1:end+n-1) = shps(six).X(1:end-1);
    rfhd.lat(end+1:end+n-1) = shps(six).Y(1:end-1);
    rfhd.polylon(end+1:end+n) = shps(six).X;
    rfhd.polylat(end+1:end+n) = shps(six).Y;

    if ( strcmp(shps(six).DESCRIPT,'Coral Reef') )
      reef.lon(end+1:end+n-1) = shps(six).X(1:end-1);
      reef.lat(end+1:end+n-1) = shps(six).Y(1:end-1);
      reef.polylon(end+1:end+n) = shps(six).X;
      reef.polylat(end+1:end+n) = shps(six).Y;
    end;
  end;

  disp(['Saving ',matfname]);
  save(matfname,'shps','rfhd','reef');
end;

clear ans coastpath matfname n six;
