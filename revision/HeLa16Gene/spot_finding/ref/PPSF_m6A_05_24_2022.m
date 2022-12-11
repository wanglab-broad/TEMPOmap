%% m6A minMax_morphRecon_SFr2_111

addpath('/stanley/WangLab/Documents/starmap_colab/Pipeline/');
addpath('/stanley/WangLab/Documents/starmap_colab/Code/matlab/');
addpath('/stanley/WangLab/Documents/starmap_colab/Code/matlab/myfunction/');

f_paths = {'/stanley/WangLab/Data/Processed/2022-04-14-m6A-100-gene/'};
run_id = '5-26-2022';
useGPU = false;

f_dims = {[2048 2048 51 4 5]}; % whole image -> make sure last 3 digits identical (dims & gblock)
f_gblock = {[1024 1024 51 4 5]};

datasets = struct('input_path', f_paths, 'input_dim', f_dims, 'gpu_block', f_gblock);
input_path = datasets.input_path;
input_dim = datasets.input_dim;
gpu_block = datasets.gpu_block;

disp('check check')

curr_data_dir = tile
curr_out_path = fullfile(input_path, 'output', run_id, curr_data_dir)
if ~exist(curr_out_path, 'dir')
    mkdir(curr_out_path)
end

%%
starting = tic;
sdata = new_STARMapDataset(input_path, 'useGPU', useGPU);
sdata.log = fopen(fullfile(curr_out_path, 'log.txt'), 'w');
fprintf(sdata.log, sprintf("====Current Tile: %s====\n", curr_data_dir));
sdata = sdata.LoadRawImages('sub_dir', curr_data_dir, 'zrange', [1:45], 'output_class', "cell");
sdata.rawImages = im_cell2mat(sdata.rawImages)
sdata = sdata.SwapChannels;
sdata = sdata.MinMaxNormalize;
sdata = sdata.MorphoRecon('Method', "2d", 'radius', 2);
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

%%
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
    bdata = bdata.SpotFinding('Method', "max3d", 'ref_index', 2, 'showPlots', false);
    bdata = bdata.ReadsExtraction('voxelSize', [1 1 1]);
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
% upd(p);

sdata.jobFinished


toc(starting) / 3600
