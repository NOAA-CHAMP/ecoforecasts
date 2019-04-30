1;
%% SCRIPT to run simple tests of the grid logic in INTERP_FIELD

fld = ones(7,8); lats=1:size(fld,1); lons=1:size(fld,2); [LON,LAT]=meshgrid(lons,lats); fld,
newfld=interp_field(lats,lons,fld,LAT,LON,{@nansum,2,2},nan); newfld=reshape(newfld,size(LAT)),

fld = ones(7,8); lats=1:size(fld,1); lons=1:size(fld,2); [LON,LAT]=meshgrid(lons,lats); fld,
newfld=interp_field(lats,lons,fld,LAT,LON,{@nansum,2,2,5},nan); newfld=reshape(newfld,size(LAT)),

fld = ones(7,8); lats=1:size(fld,1); lons=1:size(fld,2); [LON,LAT]=meshgrid(lons,lats); fld,
newfld=interp_field(lats,lons,fld,LAT,LON,{@nansum,2,2,15},nan); newfld=reshape(newfld,size(LAT)),

fld = ones(7,8); lats=1:size(fld,1); lons=1:size(fld,2); [LON,LAT]=meshgrid(lons,lats); fld,
newfld=interp_field(lats,lons,fld,LAT,LON,{@nansum,2,2,17},nan); newfld=reshape(newfld,size(LAT)),


fld = ones(7,8); lats=1:size(fld,1); lons=1:size(fld,2); [LON,LAT]=meshgrid(lons-0.5,lats-0.5); fld,
newfld=interp_field(lats,lons,fld,LAT,LON,{@nansum,2,2},nan); newfld=reshape(newfld,size(LAT)),

fld = ones(7,8); lats=1:size(fld,1); lons=1:size(fld,2); [LON,LAT]=meshgrid(lons-0.5,lats+0.5); fld,
newfld=interp_field(lats,lons,fld,LAT,LON,{@nansum,2,2},nan); newfld=reshape(newfld,size(LAT)),

fld = ones(7,8); lats=1:size(fld,1); lons=1:size(fld,2); [LON,LAT]=meshgrid(lons,lats-0.5); fld,
newfld=interp_field(lats,lons,fld,LAT,LON,{@nansum,2,2},nan); newfld=reshape(newfld,size(LAT)),

fld = ones(7,8); lats=1:size(fld,1); lons=1:size(fld,2); [LON,LAT]=meshgrid(lons+4,lats-4); fld,
newfld=interp_field(lats,lons,fld,LAT,LON,{@nansum,2,2},nan); newfld=reshape(newfld,size(LAT)),

fld = ones(7,8); lats=1:size(fld,1); lons=1:size(fld,2); [LON,LAT]=meshgrid(lons+15,lats-15); fld,
newfld=interp_field(lats,lons,fld,LAT,LON,{@nansum,2,2},nan); newfld=reshape(newfld,size(LAT)),
