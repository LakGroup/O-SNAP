save_dir = uigetdir('H:\','Select analysis directory');

if ~exist(fullfile(save_dir,"sim_data1"),'dir')
    mkdir(fullfile(save_dir,"sim_data1"))
end
if ~exist(fullfile(save_dir,"sim_data2"),'dir')
    mkdir(fullfile(save_dir,"sim_data2"))
end

generate_nucleus_set(fullfile(save_dir,"sim_data1"),80);
batch_generate_locs(data, fullfile(save_dir,"sim_data1"));

generate_nucleus_set(fullfile(save_dir,"sim_data2"),80,"major_axis",[2e3,2.4e2],"minor_axis",[1.3e3,2e2]);
batch_generate_locs(data, fullfile(save_dir,"sim_data2"));

