clear, clc, close all

project_path = "/Users/joshuadimasaka/Desktop/PhD/GitHub/riskaudit/";
path_meteor_groundtruth = project_path + ...
    "data/groundtruth/METEOR_PROJECT_2002/";

country_path = struct2table(dir(path_meteor_groundtruth));
country_path = string(table2cell(country_path(4:end,1)));

% select country 'i'
for i = 1:numel(country_path)
% for i = 1:numel(custom_list) %numel(country_path)
    tic
    disp(country_path(i))

    % set master path for the selected country
    path_i = path_meteor_groundtruth + country_path(i) + "/";
    
    % load raster grid of location IDs
    loc_raster_file = dir(fullfile(path_i, '*rectangles_rasterized.tif')).name;
    loc_raster_path = path_i + string(loc_raster_file);
    [loc,locR] = readgeoraster(loc_raster_path);

    % make new directory for generated tiles wherein the criterion is to
    % have a balanced representation of classes of vulnerability (including
    % unlabeled) so that our model can learn to distinguish fairly
    [~,~] = mkdir(path_i+"whole/extents");

    % export the extent information as geojson (for use of 
    % Google Earth/S1/S2/L8/L9 datasets)
    template_geojson = ...
        ['{ "type": "FeatureCollection", "name": ' ...
        '"XXX_country", ' ...                   % country code
        '"crs": { "type": "name", "properties": ' ...
        '{ "name": "urn:ogc:def:crs:OGC:1.3:CRS84" } }, ' ...
        '"features": [{ "type": "Feature", "properties": { }, ' ...
        '"geometry": { "type": "Polygon", ' ...
        '"coordinates": [ [ ' ...
        '[ ll_lot, ll_lat ], ' ...              % lower left
        '[ lr_lot, lr_lat ], ' ...              % lower right
        '[ ur_lot, ur_lat ], ' ...              % upper right
        '[ ul_lot, ul_lat ], ' ...              % upper left
        '[ ll_lot, ll_lat ] ] ] } }]}'];        % lower left
    ll_lot_list = [];
    ll_lat_list = [];
    lr_lot_list = [];
    lr_lat_list = [];
    ur_lot_list = [];
    ur_lat_list = [];
    ul_lot_list = [];
    ul_lat_list = [];

    ll_lot = locR.LongitudeLimits(1);
    ll_lat = locR.LatitudeLimits(1);% + locR.CellExtentInLatitude/2;
    lr_lot = locR.LongitudeLimits(2);
    lr_lat = ll_lat;
    ur_lot = lr_lot;
    ur_lat = locR.LatitudeLimits(2);% - locR.CellExtentInLatitude/2;
    ul_lot = ll_lot;
    ul_lat = ur_lat;
    ll_lot_list = [ll_lot_list; ll_lot];
    ll_lat_list = [ll_lat_list; ll_lat];
    lr_lot_list = [lr_lot_list; lr_lot];
    lr_lat_list = [lr_lat_list; lr_lat];
    ur_lot_list = [ur_lot_list; ur_lot];
    ur_lat_list = [ur_lat_list; ur_lat];
    ul_lot_list = [ul_lot_list; ul_lot];
    ul_lat_list = [ul_lat_list; ul_lat];
    writetable(table(ll_lot_list, ll_lat_list, lr_lot_list, lr_lat_list, ...
                     ur_lot_list, ur_lat_list, ul_lot_list, ul_lat_list), ...
                     path_i+loc_raster_file(1:3)+"_tile_geoextents.csv",'Delimiter',',')

    newStr = strrep(template_geojson,'XXX',loc_raster_file(1:3));
    newStr = strrep(newStr,'ll_lot', num2str(ll_lot,'%.15f'));
    newStr = strrep(newStr,'ll_lat', num2str(ll_lat,'%.15f'));
    newStr = strrep(newStr,'lr_lot', num2str(lr_lot,'%.15f'));
    newStr = strrep(newStr,'lr_lat', num2str(lr_lat,'%.15f'));
    newStr = strrep(newStr,'ur_lot', num2str(ur_lot,'%.15f'));
    newStr = strrep(newStr,'ur_lat', num2str(ur_lat,'%.15f'));
    newStr = strrep(newStr,'ul_lot', num2str(ul_lot,'%.15f'));
    newStr = strrep(newStr,'ul_lat', num2str(ul_lat,'%.15f'));

    fileID = fopen(path_i+"whole/extents/"+ ...
        loc_raster_file(1:3)+"_all"+'.geojson','w');
    fprintf(fileID, '%s\n', string(newStr)); 
    fclose(fileID);

    toc
end