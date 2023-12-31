clear, clc, close all

project_path = "/Users/joshuadimasaka/Desktop/PhD/GitHub/riskaudit/";
path_meteor_groundtruth = project_path + ...
    "data/groundtruth/METEOR_PROJECT_2002/";

country_path = struct2table(dir(path_meteor_groundtruth));
country_path = string(table2cell(country_path(4:end,1)));

% select country 'i'
for i = 1:numel(country_path)
    tic, disp(i)

    % set master path for the selected country
    path_i = path_meteor_groundtruth + country_path(i) + "/";
    
    % load raster grid of location IDs
    loc_raster_file = dir(fullfile(path_i, '*rectangles_rasterized.tif')).name;
    loc_raster_path = path_i + string(loc_raster_file);
    [loc,locR] = readgeoraster(loc_raster_path);
    
    % load csv info of floor area and building counts
    csv_file = dir(fullfile(path_i, '*loc*')).name;
    csv_path = path_i + string(csv_file);
    labels = readtable(csv_path);
    bldgtype = unique(labels.OrgConstructionCode);

    % make new directory for rasterized whole map of attributes
    [~,~] = mmkdir(path_i+"attr_rasterized");
    
    % select building type 'j'
    for j = 1:numel(bldgtype)
    
        % get indices to assign values
        loc_j = labels.LocGroup(string(labels.OrgConstructionCode)==string(bldgtype{j,1}));
        [~,Locb] = ismember(loc_j,loc(:)); %Locb contains the 'loc' indices
        [~,Loc_labels] = ismember( ...
            table(loc_j, repelem(string(bldgtype{j,1}),numel(Locb))', ...
            'VariableNames',{'Var1','Var2'}),...
            table(labels.LocGroup, string(labels.OrgConstructionCode), ...
            'VariableNames',{'Var1','Var2'}) );
        
        % create arrays
        map_farea = 0.*loc;
        map_farea(Locb) = labels.FloorArea(Loc_labels);
        map_nbldg = 0.*loc;
        map_nbldg(Locb) = labels.NumberOfBuildings(Loc_labels);
        map_bldgtiv = 0.*loc;
        map_bldgtiv(Locb) = labels.BuildingTIV(Loc_labels);
        
        % save rasters
        geotiffwrite(path_i+"attr_rasterized"+"/"+ ...
            loc_raster_file(1:3)+"_farea_"+string(bldgtype{j,1})+".tif", ...
            map_farea,locR)
        geotiffwrite(path_i+"attr_rasterized"+"/"+ ...
            loc_raster_file(1:3)+"_nbldg_"+string(bldgtype{j,1})+".tif", ...
            map_nbldg,locR)
        geotiffwrite(path_i+"attr_rasterized"+"/"+ ...
            loc_raster_file(1:3)+"_bldgtiv_"+string(bldgtype{j,1})+".tif", ...
            map_bldgtiv,locR)
    end

    toc
end