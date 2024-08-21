clear all
close all

% if czis are saved in a folder (leave empty otherwise):
path_to_folder = {};

% if multiple images are stored within 1 czi file (leave empty otherwise):
path_to_czi = {'Pathtofiles'};
im_folder = ('Pathtofolder');
save_MIPs = 1;
n_channels = 4;


%% Import CZI file and save each series


ims_405 = {};
ims_488 = {};
ims_555 = {};
if n_channels == 4
    ims_647 = {};
end

show_ims = 0;
% for multiple images stored within a single .czi file

% find image files within folder
im_files = dir(fullfile(im_folder, '*.czi'));
im_paths = fullfile(im_folder, {im_files.name});
    
    data = bfopen(path_to_czi{1});
    
    % images within data are stored in this format:
    % each row corresponds to 1 series (containing all channels)
    % within each element of the first column (1,1 2,1 3,1 etc) there is an image file that corresponds to each individual image within the stack
    % images are stored as blue image (z-stack position 1), green image (z-stack position 1), red image (z-stack position 1), blue image (z-stack position 2) etc
    
    % determine number of images within the series
    [n_series,~] = size(data);

    for mm = 1:(n_series)
    path_to_czi_1 = im_paths;
    
    %%%% import ROI file
    path_to_ROI = append(path_to_czi, '_', string(mm), '_RoiSet.zip');
    ROI_file = ReadImageJROI(path_to_ROI);
    
    %rois_poly = zeros(1024,1024);
    for kk = 1:length(ROI_file)
        roi = ROI_file{kk};
        roi_coords = roi.mnCoordinates; % coordinates of the individual ROI
        im_roi = poly2mask(roi_coords(:,1), roi_coords(:,2), 1024, 1024); % create mask of each polygon
        rois_555{kk,mm} = im_roi;
    end
    %rois_488{mm} = rois_poly;
    %%%% end import ROI file

    end
    
    
    for jj = 1:n_series
        
        series = data{jj,1};
        series_names = series(:,2);
    
        
        % find the images corresponding to each channel
        if n_channels == 3
            idx_405 = find(~cellfun(@isempty,(strfind(series_names, 'C=1/3'))));
            idx_488 = find(~cellfun(@isempty,(strfind(series_names, 'C=2/3'))));
            idx_555 = find(~cellfun(@isempty,(strfind(series_names, 'C=3/3'))));
        elseif n_channels == 4
            idx_405 = find(~cellfun(@isempty,(strfind(series_names, 'C=1/4'))));
            idx_488 = find(~cellfun(@isempty,(strfind(series_names, 'C=2/4'))));
            idx_555 = find(~cellfun(@isempty,(strfind(series_names, 'C=3/4'))));
            idx_647 = find(~cellfun(@isempty,(strfind(series_names, 'C=4/4'))));
        end
        
        % move the images for each channel into their own stack
        im_stack_405 = series(idx_405,1);
        im_stack_488 = series(idx_488,1);
        im_stack_555 = series(idx_555,1);
        
        % convert from cell to 3D matrix
        im_stack_405 = cat(3,im_stack_405{:});
        im_stack_488 = cat(3,im_stack_488{:});
        im_stack_555 = cat(3,im_stack_555{:});
        
        % make MIP
        mip_405 = max(im_stack_405, [], 3);
        mip_488 = max(im_stack_488, [], 3);
        mip_555 = max(im_stack_555, [], 3);
        
        if n_channels == 4
            im_stack_647 = series(idx_647,1);
            im_stack_647 = cat(3,im_stack_647{:});
            mip_647 = max(im_stack_647, [], 3);
        end
        
        % store images in MATLAB
        index = length(ims_405) + 1;
        ims_405{index} = mip_405;
        ims_488{index} = mip_488;
        ims_555{index} = mip_555;
        if n_channels == 4
            ims_647{index} = mip_647;
        end

         %extract metadata and create metadata object for save
             metadata_im = createMinimalOMEXMLMetadata(mip_405);
             metadata = data{1,4};
             pixelSizeX = metadata.getPixelsPhysicalSizeX(0).value();
             pixelSizeX = pixelSizeX.doubleValue(); % size of each pixel (X dimension) in microns
             pixelSizeY = metadata.getPixelsPhysicalSizeY(0).value();
             pixelSizeY = pixelSizeY.doubleValue(); % size of each pixel (X dimension) in microns
             area_conv{index} = pixelSizeX * pixelSizeY; % conversion from pixels to microns (microns^2/pixels^2)
             
             pixelSize = ome.units.quantity.Length(java.lang.Double(pixelSizeX), ome.units.UNITS.MICROMETER);
             metadata_im.setPixelsPhysicalSizeX(pixelSize, 0);
             metadata_im.setPixelsPhysicalSizeY(pixelSize, 0);
             
            [px_x, px_y] = size(mip_405);
             length_x = px_x * pixelSizeX;
             length_y = px_y * pixelSizeY;
        
        % write MIP to folder
        experiment_name = extractAfter(path_to_czi{1}, 'DHRRS3_v5_noOA/'); % switch from airy to 880
        experiment_name = extractBefore(experiment_name, '.czi');
         group_index{index,1} = experiment_name;
        
        folder_path = extractBefore(path_to_czi{1}, experiment_name);
        cd(folder_path)
        
        if ~isfolder(experiment_name)
            mkdir(experiment_name)
        end

        cd(experiment_name)
        im_name = append('max_', experiment_name, '_n', string(jj));
        
        % write MIPs to folder
        imwrite(mip_405, append(im_name, '_405.tiff'))
        imwrite(mip_488, append(im_name, '_488.tiff'))
        imwrite(mip_555, append(im_name, '_555.tiff'))
        if n_channels == 4
            imwrite(mip_647, append(im_name, '_647.tiff'));
        end
    end
%% Background subtraction with mean as background, store files

n_images = length(ims_405);

ims_bsub_405 = {};
ims_bsub_488 = {};
ims_bsub_555 = {};
if n_channels == 4
    ims_bsub_647 = {};
end


for kk = 1:n_images
    
    % 405
    im_405 = ims_405{kk};
    mean_405 = mean(im_405(:));
    ims_bsub_405{kk} = im_405 - mean_405; %background subtraction
    
    % 488
    im_488 = ims_488{kk};
    mean_488 = mean(im_488(:));
    ims_bsub_488{kk} = im_488 - mean_488; %background subtraction
    
    
    % 555
    im_555 = ims_555{kk};
    mean_555 = mean(im_555(:));
    ims_bsub_555{kk} = im_555 - mean_555; % background subtraction
    
    
    if n_channels == 4
        im_647 = ims_647{kk};
        mean_647 = mean(im_647(:));
        ims_bsub_647{kk} = im_647 - mean_647; %background subtraction
    end
    
end

%% Threshold the LDs

n_std_for_thresh_647 = 5; %number of standard deviations to use for threshold for LipidTox

for dd = 1:n_images
    
    %segment the LDs
    im_647 = ims_bsub_647{dd};
    mean_647 = mean(im_647(:)); 
    mean_647_list{dd} = mean_647; 
    std_647 = std(double(im_647(:))); 
    figure; imshow(im_647, []); 
    
    %apply the threshold to the LDs
    thresh_647 = mean_647 + (std_647*n_std_for_thresh_647); %set adaptive threshold
    im_647_thresh{dd} = im_647 > thresh_647; %threshold the LDs
    im_647_thresh{dd} = bwareafilt(im_647_thresh{dd}, [5, 100000]);
    figure; imshow(im_647_thresh{dd}); 
    
    
end

%% Quantify LD area and intensity per cell


nuclei_index = [];

for ii = 1:n_images

    n_cells = length(find(~cellfun(@isempty, rois_555(:,ii))));
    
    for jj = 1:n_cells
        
        cell_roi = rois_555{jj,ii};

%         cell_im = false(size(ims_405{ii})); % make a blank image
%         cell_im(cyto_cc.PixelIdxList{jj}) = true;

        % segment lipid droplets
        LD_im{jj, ii} = cell_roi .* im_647_thresh{ii};
        LD_im_save = figure; imshow(LD_im{jj, ii}); 
        saveas(LD_im_save, append(experiment_name, num2str(ii), '_', num2str(jj), '_LD_thresh.tif'));
          
        %calculate the total LD area per cell
        LD_area_pix = cellfun(@nnz, LD_im); 
        LD_area_pix = num2cell(LD_area_pix); 
        LD_area_micron{jj, ii} = LD_area_pix{jj, ii} * area_conv{ii};
        writecell(LD_area_micron, 'LD_area_micron.xls'); 
        
        %cacluate the average LD intensity per cell
        LD_ints = uint16(double(ims_bsub_647{ii}) .* double(LD_im{jj, ii})); 
        LD_ints = double(LD_ints);
        LD_ints(LD_ints == 0) = NaN;
        LD_ints_mean{jj,ii} = mean(LD_ints(:), 'omitnan'); 
        writecell(LD_ints_mean, 'LD_ints_mean.xls'); 
    end
end


close all


