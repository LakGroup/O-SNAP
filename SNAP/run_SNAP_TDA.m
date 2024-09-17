%% python setup
py_env_version = '3.11';
py_env_dir="C:\Users\HannahKim\.pyenv\pyenv-win\versions\3.11.9\python.exe";
if exist('py_env', 'var')
    if py_env.Status == "Loaded" && strcmp(py_env.Version, py_env_version)
        disp(['To change the Python version, restart MATLAB, then call pyenv(Version="' py_env_version '").'])
    end
else
    py_env = pyenv(Version=py_env_dir);
end

work_dir = fullfile(root_dir,analysis_name);

%% load data
data_info_table = get_valid_voronoi_data(work_dir,groups,reps,{'x','y'});
filepaths = data_info_table{:,'filepath'}; % TODO: MAKE RUN PARALLEL
filepaths = filepaths(1:2); % TODO: DELETE

%% run code
output_data = cell(size(filepaths,1),1);
for i=1:size(filepaths,1)
    data = load_variables(filepaths(i),{'x','y'});
    output_data{i} = double(pyrunfile("C:\Users\HannahKim\Documents\STORM\SNAP\run_SNAP_TDA.py","output_data",x=data.x,y=data.y));
end
% TODO: remove where output_data = 0

%% plot
output = output_data{1};
birth = output(:,1);
death = output(:,2);

figure();
tiledlayout('flow')
nexttile
scatter(birth,death);
xlabel("birth [nm]")
ylabel("death [nm]")
nexttile
for i=1:size(output,1)
    plot([death death-birth],[size(output,1)-i size(output,1)-i],'b')
    hold on
end
xlabel("birth [nm]")
ylabel("hole #")