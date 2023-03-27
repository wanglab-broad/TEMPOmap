tile='Position015'
%% TEMPOmap 16Gene

addpath('/stanley/WangLab/Documents/starmap_colab/Pipeline/');
addpath('/stanley/WangLab/Documents/starmap_colab/Code/matlab/');
addpath('/stanley/WangLab/Documents/starmap_colab/Code/matlab/myfunction/');

input_path = '/stanley/WangLab/Data/Processed/2022-09-12-Rena-HeLa16Gene/';
run_id = '2022-09-14'

input_dim = [2048 2048 35 4 3];
gpu_block = [1024 1024 35 4 3];

useGPU = false;

curr_data_dir = tile
curr_out_path = fullfile(input_path, 'output', run_id, curr_data_dir)

if ~exist(curr_out_path, 'dir')
    mkdir(curr_out_path)
end

starting = tic;
sdata = new_STARMapDataset(input_path, 'useGPU', useGPU);
sdata.log = fopen(fullfile(curr_out_path, 'log.txt'), 'w');
fprintf(sdata.log, sprintf("====Current Tile: %s====\n", curr_data_dir));
sdata = sdata.LoadRawImages('sub_dir', curr_data_dir, 'input_dim', input_dim);
sdata = sdata.SwapChannels; % !!
%% TO ADD: min-max (ref. immune res)

sdata = sdata.HistEqualize('Method', "inter_round");
sdata = sdata.HistEqualize('Method', "intra_round");
sdata = sdata.MorphoRecon('Method', "2d", 'radius', 3); % param-to-tune: radius

sdata = sdata.test_GlobalRegistration;

gpu_overlap = 0;
gb_idx = BlockTest(sdata.rawImages, gpu_block, gpu_overlap);
NgpuBlock = numel(gb_idx)

tile_goodSpots = [];
tile_goodReads = [];

tile_allSpots = [];
tile_allReads = [];

tile_score = [];
tile_counts = 0;

for i=1:NgpuBlock

    bdata_start = tic;
    gpu_idx = gb_idx{i}

    fprintf(sdata.log, sprintf("====Current Block: %d====\n", i));
    fprintf(sdata.log, sprintf("Index: %d - %d - %d - %d\n", gpu_idx));

    bdata = new_STARMapDataset(input_path, 'useGPU', useGPU);
    bdata.log = sdata.log;
    bdata = bdata.LoadDim(gpu_block);
    bdata.registeredImages = sdata.registeredImages(...
        gpu_idx(1,1):gpu_idx(1,2),...
        gpu_idx(2,1):gpu_idx(2,2),...
        :,:,:);

    bdata = bdata.LocalRegistration('Iterations', 50, 'AccumulatedFieldSmoothing', 1);

    % replace the old image
    sdata.registeredImages(...
        gpu_idx(1,1):gpu_idx(1,2),...
        gpu_idx(2,1):gpu_idx(2,2),...
        :,:,:) = bdata.registeredImages;

    bdata = bdata.LoadCodebook;
    bdata = bdata.SpotFinding('Method', "max3d", 'ref_index', 1, 'showPlots', false);
    bdata = bdata.ReadsExtraction('voxelSize', [2 2 1]);
    bdata = bdata.ReadsFiltration('showPlots', false);

    % block_offsets = [gpu_idx(2,1) - 1 gpu_idx(1,1) - 1 0]
    if size(bdata.goodSpots, 1) ~= 0
        block_offsets = [gpu_idx(2,1) - 1 gpu_idx(1,1) - 1 0]
        block_offsets = repmat(block_offsets, size(bdata.goodSpots, 1), 1);
        bdata.goodSpots = bdata.goodSpots + int16(block_offsets);

        block_offsets = [gpu_idx(2,1) - 1 gpu_idx(1,1) - 1 0]
        block_offsets = repmat(block_offsets, size(bdata.allSpots, 1), 1);
        bdata.allSpots = bdata.allSpots + int16(block_offsets);
    end

    % Update tile reads
    tile_goodSpots = [tile_goodSpots; bdata.goodSpots];
    tile_goodReads = [tile_goodReads; bdata.goodReads];

    tile_allSpots = [tile_allSpots; bdata.allSpots];
    tile_allReads = [tile_allReads; bdata.allReads];

    tile_score = [tile_score; bdata.FilterScores];
    tile_counts = tile_counts + size(bdata.allSpots, 1);

    fprintf(sprintf("====Block %d Finished [time=%02f]====", i, toc(bdata_start)));
end

% Save round 1 image
output_dir = fullfile(input_path, 'output', run_id, 'round1_merged');
if ~exist(output_dir, 'dir')
   mkdir(output_dir);
end
r1_img = max(sdata.registeredImages(:,:,:,:,1), [], 4);
r1_img_name = fullfile(output_dir, sprintf("%s.tif", curr_data_dir));
SaveSingleTiff(r1_img, r1_img_name)

% Save tile points
save(fullfile(curr_out_path, strcat('goodPoints_max3d.mat')), 'tile_goodReads', 'tile_goodSpots');
save(fullfile(curr_out_path, strcat('allPoints_max3d.mat')), 'tile_allReads', 'tile_allSpots');

fprintf(sdata.log, sprintf("Average Score: %.2f - %.2f - %.2f\n", mean(tile_score, 'omitnan')));

% Save dots image
curr_img = max(sdata.registeredImages(:,:,:,:,1), [], 4);
curr_img = max(curr_img, [], 3);

curr_img_out_path = fullfile(input_path, 'output', run_id, 'dots_image')

if ~exist(curr_img_out_path, 'dir')
    mkdir(curr_img_out_path)
end

if ~isempty(tile_goodSpots)
    plot_centroids(tile_goodSpots, curr_img, 2, 'r')
    saveas(gcf, fullfile(curr_img_out_path, sprintf("%s.tif", curr_data_dir)));
end

reset(gpuDevice)
fclose(sdata.log);
sdata.jobFinished

toc(starting) / 3600
