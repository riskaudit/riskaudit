clear, clc, close all

project_path = "/Users/joshuadimasaka/Desktop/PhD/GitHub/riskaudit/";
path_meteor_groundtruth = project_path + ...
    "data/groundtruth/METEOR_PROJECT_2002/";

country_path = struct2table(dir(path_meteor_groundtruth));
country_path = string(table2cell(country_path(4:end,1)));

tile_edge = 8; % tile_edge = 8 is approx. 4km by 4km area
ntiles_per_country = 100; %random sampling - not practical to use all of the map (size issues) - must be multiples of 10
seed = 1; %reproducibility

% select country 'i'
% for i = 1:numel(country_path)
for i = 1:numel(country_path)
    tic
    disp(country_path(i))

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

    % make new directory for generated tiles wherein the criterion is to
    % have a balanced representation of classes of vulnerability (including
    % unlabeled) so that our model can learn to distinguish fairly
    [~,~] = mkdir(path_i+"tiles/extents");
    [~,~] = mkdir(path_i+"tiles/images");
    
    % initialize map_nlabel
    map_nlabel = zeros(size(loc,1),size(loc,2),numel(bldgtype));

    % select building type 'j'
    for j = 1:numel(bldgtype)
        % map that shows how many labels are present
        map_nlabel(:,:,j) = double(readgeoraster(path_i+...
            "attr_rasterized"+"/"+ ...
            loc_raster_file(1:3)+"_farea_"+...
            string(bldgtype{j,1})+".tif")>0);
    end

    % load country mask
    countryMask_file = dir(fullfile(path_i, '*_country.tif')).name;
    countryMask_path = path_i + string(countryMask_file);
    [countryMask,~] = readgeoraster(countryMask_path);
    valid_col = find(sum(countryMask,1)~=0)';
    valid_row = find(sum(countryMask,2)~=0);
    valid_cole = 1:(size(loc,2)-tile_edge+1);
    valid_rowe = 1:(size(loc,1)-tile_edge+1);
    valid_colF = valid_cole(ismember(valid_cole,valid_col))';
    valid_rowF = valid_rowe(ismember(valid_rowe,valid_row))';
    

    % determine number of tiles and its geographic extent (having balanced
    % representation of all labels) - greedy dynamic search
    ind = [];
    attr = [];
    nzeros_ratio = [];
    prop = 0.5; % proportion of zeros
    cover = 1;
    overlap_ratio = 0;
    tic
    for icol = 1:numel(valid_colF)
        icol*100/numel(valid_colF)
        col = valid_colF(icol);
        for irow = 1:numel(valid_rowF)
            row = valid_rowF(irow);
            if sum(sum(countryMask(row:(row+tile_edge-1),col:(col+tile_edge-1))>0)) ~= 0
                % cropped to selection
                tile = map_nlabel(row:(row+tile_edge-1),col:(col+tile_edge-1),:);
                % linear indices
                ind_rc = sub2ind(size(loc),...
                       repmat(row:(row+tile_edge-1),[1 tile_edge]),...
                       repelem(col:(col+tile_edge-1),tile_edge));
                % N - number of classes present (integer)
                N = numel(unique(reshape(sum(tile,3),1,[])'));
                % nzeros - number of zeros or unlabeled segments
                nzeros_temp = sum(reshape(sum(tile,3),1,[])' == 0);
                % multi-class criterion
                attr_temp = reshape(sum(sum(tile,1),2),1,[]);
                if N > 1 && ... % ensures 2 or more classes
                   sum(ismember(ind_rc(:),ind(:))) == overlap_ratio*(tile_edge^2) && ... % ensures no overlapping of tiles
                   sum(sum(countryMask(row:(row+tile_edge-1),col:(col+tile_edge-1))>0)) == cover*(tile_edge^2) && ... % ensures within country borders
                   nzeros_temp < prop*(tile_edge^2) % ensures (1) computational efficiency & 
                    % (2) we would like to train our machine learning model to 
                    % learn the pattern half from those with labels and half 
                    % from unlabeled segments. We do not need those with, 
                    % let's say, 75% or 90% unlabeled segments, because the 
                    % future steps will need to incorporate the prior belief 
                    % from meteor maps that those areas have no labels unless 
                    % significant patterns start to show from auxialliary 
                    % datasets like Landsat or Sentinel
                    ind = [ind; ind_rc];
                    attr = [attr; attr_temp];
                    nzeros_ratio = [nzeros_ratio; nzeros_temp/(tile_edge^2)];
                end
            end

        end
    end
    toc

    
    % 'importance sampling'
    if size(ind,1) ~= 0
        importance_weight = 0.*attr;
        for y = 1:numel(bldgtype)
            inter_x = histogram(attr(:,y),[0 1:1:tile_edge]*tile_edge,'Normalization','pdf').BinEdges;
            inter_pdf = histogram(attr(:,y),[0 1:1:tile_edge]*tile_edge,'Normalization','pdf').Values;
            for z = 1:size(attr,1)
                if attr(z,y) == 0
                    importance_weight(z,y) = inter_pdf(1);
                else
                    importance_weight(z,y) = inter_pdf(find((attr(z,y) <= inter_x)==1,1,'first')-1);
                end
            end
        end
        % convert into single weight (assuming independence)
        importance_weight_single = ones(size(attr,1),1);
        for y = 1:numel(bldgtype)
            importance_weight_single = importance_weight_single.*importance_weight(:,y);
        end
        % select the tiles to have balanced representation - our criteria is
        % the number of zeros (no labels) in a given tile. Fewer zeros mean
        % that those tiles are from highly condensed areas whereas tiles with
        % many zeros are from areas that could possible be remote or small
        % communities.
        selected_rows = [];
        n_selected = [];
        k_selected = [];
        for w = 1:10
            w
            idx = find(round((nzeros_ratio*(1/prop)*10+prop))==w);
            n_selected = [n_selected; numel(idx)];
        end
    end

    % second fix - number of extracted tiles is enough but selected rows
    % are just limiting
    if sum(n_selected.*(n_selected<=10)) < ntiles_per_country && size(ind,1) > ntiles_per_country
        s = find((n_selected.*(n_selected<=10))==0 & n_selected~=0);
        [B,I] = sort(n_selected(s));
        running_value = sum(n_selected.*(n_selected<=10));
        distribution_value = round((ntiles_per_country-running_value)./numel(s));
        for x = 1:numel(s)
            if n_selected(s(I(x))) < distribution_value
                n_selected(s(I(x))) = n_selected(s(I(x)));
                running_value = running_value + n_selected(s(I(x)));
                distribution_value = round((ntiles_per_country-running_value)./(numel(s)-x));
            elseif n_selected(s(I(x))) >= distribution_value
                running_value = running_value + distribution_value;
                n_selected(s(I(x))) = distribution_value;
            end
        end
        if sum(n_selected) ~= ntiles_per_country 
           n_selected(end) = n_selected(end) - (sum(n_selected)-ntiles_per_country);
        end
        for w = 1:10
            idx = find(round((nzeros_ratio*(1/prop)*10+prop))==w);
            [~,I] = maxk(importance_weight_single(idx),n_selected(w));
            selected_rows = [selected_rows; idx(I)];
        end
        ind = ind(selected_rows,:);
    elseif size(ind,1) < ntiles_per_country
        % third fix - number of extracted files is not enough ()
        while size(ind,1) < ntiles_per_country && cover >= 0.0001
            tic
            prop = 0.5;
            while size(ind,1) < ntiles_per_country && prop < 1
                for icol = 1:numel(valid_colF)
                    icol*100/numel(valid_colF)
                    col = valid_colF(icol);
                    for irow = 1:numel(valid_rowF)
                        row = valid_rowF(irow);
                        if sum(sum(countryMask(row:(row+tile_edge-1),col:(col+tile_edge-1))>0)) ~= 0
                            tile = map_nlabel(row:(row+tile_edge-1),col:(col+tile_edge-1),:);
                            ind_rc = sub2ind(size(loc),...
                                   repmat(row:(row+tile_edge-1),[1 tile_edge]),...
                                   repelem(col:(col+tile_edge-1),tile_edge));
                            N = numel(unique(reshape(sum(tile,3),1,[])'));
                            nzeros_temp = sum(reshape(sum(tile,3),1,[])' == 0);
                            attr_temp = reshape(sum(sum(tile,1),2),1,[]);
    
                            if N > 1 && ... % ensures 2 or more classes
                               sum(ismember(ind_rc(:),ind(:))) == overlap_ratio*(tile_edge^2) && ... % ensures no overlapping of tiles
                               sum(sum(countryMask(row:(row+tile_edge-1),col:(col+tile_edge-1))>0)) >= cover*(tile_edge^2) && ... % ensures within country borders
                               nzeros_temp < prop*(tile_edge^2)
    
                                ind = [ind; ind_rc];
                                attr = [attr; attr_temp];
                                nzeros_ratio = [nzeros_ratio; nzeros_temp/(tile_edge^2)];
                            end
                        end
                    end
                end
                prop = prop + 0.02; % proportion of zeros
            end
            cover = cover - 0.05
            toc
        end


        importance_weight = 0.*attr;
        for y = 1:numel(bldgtype)
            inter_x = histogram(attr(:,y),[0 1:1:tile_edge]*tile_edge,'Normalization','pdf').BinEdges;
            inter_pdf = histogram(attr(:,y),[0 1:1:tile_edge]*tile_edge,'Normalization','pdf').Values;
            for z = 1:size(attr,1)
                if attr(z,y) == 0
                    importance_weight(z,y) = inter_pdf(1);
                else
                    importance_weight(z,y) = inter_pdf(find((attr(z,y) <= inter_x)==1,1,'first')-1);
                end
            end
        end
        importance_weight_single = ones(size(attr,1),1);
        for y = 1:numel(bldgtype)
            importance_weight_single = importance_weight_single.*importance_weight(:,y);
        end
        
        selected_rows = [];
        n_selected = [];
        k_selected = [];
        for w = 1:10
            % find(round((nzeros_ratio*(1/prop)*10+prop))==w);
            idx = find(ceil((nzeros_ratio+0.001)*10)==w);
            n_selected = [n_selected; numel(idx)];
        end
        if sum(n_selected.*(n_selected<=10)) < ntiles_per_country
            s = find((n_selected.*(n_selected<=10))==0 & n_selected~=0);
            [B,I] = sort(n_selected(s));
            running_value = sum(n_selected.*(n_selected<=10));
            distribution_value = round((ntiles_per_country-running_value)./numel(s));
            for x = 1:numel(s)
                if n_selected(s(I(x))) < distribution_value
                    n_selected(s(I(x))) = n_selected(s(I(x)));
                    running_value = running_value + n_selected(s(I(x)));
                    distribution_value = round((ntiles_per_country-running_value)./(numel(s)-x));
                elseif n_selected(s(I(x))) >= distribution_value
                    running_value = running_value + distribution_value;
                    n_selected(s(I(x))) = distribution_value;
                end
            end
            if sum(n_selected) ~= ntiles_per_country 
               n_selected(end) = n_selected(end) - (sum(n_selected)-ntiles_per_country);
            end
            for w = 1:10
                idx = find(ceil((nzeros_ratio+0.001)*10)==w);
                [~,I] = maxk(importance_weight_single(idx),n_selected(w));
                selected_rows = [selected_rows; idx(I)];
            end
            ind = ind(selected_rows,:);
        end
        
    end


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
    for k = 1:size(ind,1)

        [row,col] = ind2sub(size(loc),ind(k,:));
        ll_lot = locR.LongitudeLimits(1) + (min(col)-1).*locR.CellExtentInLongitude;
        ll_lat = locR.LatitudeLimits(2) - max(row).*locR.CellExtentInLatitude;
        lr_lot = locR.LongitudeLimits(1) + max(col).*locR.CellExtentInLongitude;
        lr_lat = locR.LatitudeLimits(2) - max(row).*locR.CellExtentInLatitude;
        ur_lot = locR.LongitudeLimits(1) + max(col).*locR.CellExtentInLongitude;
        ur_lat = locR.LatitudeLimits(2) - (min(row)-1).*locR.CellExtentInLatitude;
        ul_lot = locR.LongitudeLimits(1) + (min(col)-1).*locR.CellExtentInLongitude;
        ul_lat = locR.LatitudeLimits(2) - (min(row)-1).*locR.CellExtentInLatitude;
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

        fileID = fopen(path_i+"tiles/extents/"+ ...
            loc_raster_file(1:3)+"_"+num2str(k)+"_of_"+num2str(size(ind,1))+"_"+...
            num2str(min(row))+"_"+...
            num2str(max(row))+"_"+...
            num2str(min(col))+"_"+...
            num2str(max(col))+...
            '.geojson','w');
        fprintf(fileID, '%s\n', string(newStr)); 
        fclose(fileID);

        for m = 1:numel(bldgtype)

            [map_farea_crop,locR_crop] = geocrop(...
                readgeoraster(path_i+...
                "attr_rasterized"+"/"+ ...
                loc_raster_file(1:3)+"_farea_"+...
                string(bldgtype{m,1})+".tif"),locR, ...
                [ceil(ll_lat*1000)/1000 floor(ul_lat*1000)/1000], ...
                [ceil(ll_lot*1000)/1000 floor(lr_lot*1000)/1000]);
            geotiffwrite(path_i+"tiles/images/"+ ...
                loc_raster_file(1:3)+"_farea_"+string(bldgtype{m,1})+"_"+ ...
                num2str(k)+"_of_"+num2str(size(ind,1))+"_"+...
                num2str(min(row))+"_"+...
                num2str(max(row))+"_"+...
                num2str(min(col))+"_"+...
                num2str(max(col))+".tif", ...
                map_farea_crop, ...
                locR_crop);

            [map_nbldg_crop,~] = geocrop(...
                readgeoraster(path_i+...
                "attr_rasterized"+"/"+ ...
                loc_raster_file(1:3)+"_nbldg_"+...
                string(bldgtype{m,1})+".tif"),locR, ...
                [ceil(ll_lat*1000)/1000 floor(ul_lat*1000)/1000], ...
                [ceil(ll_lot*1000)/1000 floor(lr_lot*1000)/1000]);
            geotiffwrite(path_i+"tiles/images/"+ ...
                loc_raster_file(1:3)+"_nbldg_"+string(bldgtype{m,1})+"_"+ ...
                num2str(k)+"_of_"+num2str(size(ind,1))+"_"+...
                num2str(min(row))+"_"+...
                num2str(max(row))+"_"+...
                num2str(min(col))+"_"+...
                num2str(max(col))+".tif", ...
                map_nbldg_crop, ...
                locR_crop);

            [map_bldgtiv_crop,~] = geocrop(...
                readgeoraster(path_i+...
                "attr_rasterized"+"/"+ ...
                loc_raster_file(1:3)+"_bldgtiv_"+...
                string(bldgtype{m,1})+".tif"),locR, ...
                [ceil(ll_lat*1000)/1000 floor(ul_lat*1000)/1000], ...
                [ceil(ll_lot*1000)/1000 floor(lr_lot*1000)/1000]);
            geotiffwrite(path_i+"tiles/images/"+ ...
                loc_raster_file(1:3)+"_bldgtiv_"+string(bldgtype{m,1})+"_"+ ...
                num2str(k)+"_of_"+num2str(size(ind,1))+"_"+...
                num2str(min(row))+"_"+...
                num2str(max(row))+"_"+...
                num2str(min(col))+"_"+...
                num2str(max(col))+".tif", ...
                map_bldgtiv_crop, ...
                locR_crop);

            if      size(map_farea_crop,1) ~= tile_edge
                k
                disp("cropped map not correct size")
            elseif  size(map_nbldg_crop,1) ~= tile_edge
                k
                disp("cropped map not correct size")
            elseif  size(map_bldgtiv_crop,1) ~= tile_edge
                k
                disp("cropped map not correct size")
            end
        end
    end
    toc


end