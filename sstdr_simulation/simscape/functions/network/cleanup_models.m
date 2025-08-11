function cleanup_models()
%CLEANUP_MODELS Clean up scattered Simulink model files
%
% This function finds all .slx, .slxc, and .autosave files in the simscape
% directory and its subdirectories (except models/) and either deletes them
% or moves them to the models/ directory.
%
% Usage:
%   cleanup_models()

fprintf('=== Cleaning up scattered model files ===\n');

% Get current directory
current_dir = pwd;
simscape_root = fileparts(mfilename('fullpath'));

% Change to simscape root directory
cd(simscape_root);

% Ensure models directory exists
models_dir = 'models';
if ~exist(models_dir, 'dir')
    mkdir(models_dir);
    fprintf('Created models/ directory\n');
end

% Find all model files (excluding models/ directory)
fprintf('Searching for scattered model files...\n');

% Search patterns
patterns = {'*.slx', '*.slxc', '*.autosave'};
scattered_files = {};

for i = 1:length(patterns)
    % Find files in root directory
    files = dir(patterns{i});
    for j = 1:length(files)
        if ~files(j).isdir
            scattered_files{end+1} = files(j).name;
        end
    end
    
    % Find files in subdirectories (except models/)
    subdirs = dir('.');
    for j = 1:length(subdirs)
        if subdirs(j).isdir && ~strcmp(subdirs(j).name, '.') && ...
           ~strcmp(subdirs(j).name, '..') && ~strcmp(subdirs(j).name, 'models')
            subdir_files = dir(fullfile(subdirs(j).name, patterns{i}));
            for k = 1:length(subdir_files)
                if ~subdir_files(k).isdir
                    scattered_files{end+1} = fullfile(subdirs(j).name, subdir_files(k).name);
                end
            end
        end
    end
end

% Process found files
if isempty(scattered_files)
    fprintf('✓ No scattered model files found\n');
else
    fprintf('Found %d scattered model files:\n', length(scattered_files));
    
    deleted_count = 0;
    moved_count = 0;
    
    for i = 1:length(scattered_files)
        file_path = scattered_files{i};
        [file_dir, file_name, file_ext] = fileparts(file_path);
        
        fprintf('  %s\n', file_path);
        
        % Decide what to do with the file
        if strcmp(file_ext, '.autosave')
            % Delete autosave files
            try
                delete(file_path);
                fprintf('    ✓ Deleted (autosave)\n');
                deleted_count = deleted_count + 1;
            catch ME
                fprintf('    ✗ Could not delete: %s\n', ME.message);
            end
            
        elseif strcmp(file_ext, '.slxc')
            % Delete compiled files
            try
                delete(file_path);
                fprintf('    ✓ Deleted (compiled)\n');
                deleted_count = deleted_count + 1;
            catch ME
                fprintf('    ✗ Could not delete: %s\n', ME.message);
            end
            
        elseif strcmp(file_ext, '.slx')
            % Ask user what to do with .slx files
            response = input(sprintf('    Move to models/ directory? (y/n/d=delete) [y]: '), 's');
            if isempty(response)
                response = 'y';
            end
            
            if strcmpi(response, 'y')
                % Move to models directory
                new_path = fullfile(models_dir, [file_name file_ext]);
                if exist(new_path, 'file')
                    fprintf('    ⚠ File already exists in models/: %s\n', [file_name file_ext]);
                    overwrite = input('    Overwrite? (y/n) [n]: ', 's');
                    if ~strcmpi(overwrite, 'y')
                        fprintf('    - Skipped\n');
                        continue;
                    end
                end
                
                try
                    movefile(file_path, new_path);
                    fprintf('    ✓ Moved to models/\n');
                    moved_count = moved_count + 1;
                catch ME
                    fprintf('    ✗ Could not move: %s\n', ME.message);
                end
                
            elseif strcmpi(response, 'd')
                % Delete file
                try
                    delete(file_path);
                    fprintf('    ✓ Deleted\n');
                    deleted_count = deleted_count + 1;
                catch ME
                    fprintf('    ✗ Could not delete: %s\n', ME.message);
                end
            else
                fprintf('    - Skipped\n');
            end
        end
    end
    
    fprintf('\n=== Cleanup Summary ===\n');
    fprintf('Files processed: %d\n', length(scattered_files));
    fprintf('Files deleted: %d\n', deleted_count);
    fprintf('Files moved to models/: %d\n', moved_count);
end

% Restore original directory
cd(current_dir);

fprintf('✓ Cleanup complete\n');

end 