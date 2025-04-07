function f = loadandprocessdata(direc,varargin)
% loadandprocessdata  reads the data from the given directory from the output 
% of FastCycles 2D
%
%   f = loadandprocessdata(direct) reads the velocity V of size time_step x
%   space_points
%
%   f = loadandprocessdata(direc,'St') reads the elastic shear traction of 
%   size time_step x space_points.
%
%   The different option are 
%   the velocity: 'V'
%   the state variable: 'theta'
%   the elastic shear traction: 'St'
%   the elastic normal traction: 'Sn'
%
%   f = loadandprocessdata(direc,'V','St') reads the velocity and the 
%   elastic shear traction, both of size time_step x space_points 

% Read all files name in the directory
files = dir(direc);
file_names =  {files(:).name};

% Count the number of files
index = strfind(file_names,'output');
index = find(not(cellfun('isempty', index)));
filenames = file_names(index);
nb_files = length(filenames);


% Read the time
temp = h5read([direc filenames{1}],'/time');
time_size = length(temp);
f.time = squeeze(zeros(1,time_size*nb_files));
for i=1:length(filenames)
    temp = h5read([direc filenames{i}],'/time');
    f.time((i-1)*time_size+1:i*time_size) = temp;
end

% Find size of array
time_size_full = length(f.time);


% Read values
if ~isempty(varargin)
    for j=1:length(varargin)
        switch lower(varargin{j})
            case {'v'}
                % Find size to allocate
                temp = h5read([direc filenames{i}],'/V');
                xsize = size(temp,2);

                % Allocate
                f.V = zeros(time_size_full,xsize);

                % Read V
                for i=1:length(filenames)
                    disp(['number of files loaded: ' num2str(i)])
                    temp = h5read([direc filenames{i}],'/V');

                    f.V((i-1)*time_size+1:(i)*time_size,:) = temp;
                end

            case {'p'}
                % Find size to allocate
                temp = h5read([direc filenames{i}],'/P');
                xsize = size(temp,2);

                % Allocate
                f.P = zeros(time_size_full,xsize);

                % Read V
                for i=1:length(filenames)
                    disp(['number of files loaded: ' num2str(i)])
                    temp = h5read([direc filenames{i}],'/P');

                    f.P((i-1)*time_size+1:(i)*time_size,:) = temp;
                end


            case {'theta','t','th'}
                % Find size to allocate
                temp = h5read([direc filenames{i}],'/theta');
                xsize = size(temp,2);

                % Allocate
                f.theta = zeros(time_size_full,xsize);

                % Read theta
                for i=1:length(filenames)
                    disp(['number of files loaded: ' num2str(i)])
                    temp = h5read([direc filenames{i}],'/theta');
                    f.theta((i-1)*time_size+1:(i)*time_size,:) = temp;
                end

            case {'stressn','sn'}
                % Find size to allocate
                temp = h5read([direc filenames{i}],'/tractionNel');
                xsize = size(temp,2);

                % Allocate
                f.Sn = zeros(time_size_full,xsize);

                % Read normal traction
                for i=1:length(filenames)
                    disp(['number of files loaded: ' num2str(i)])
                    temp = h5read([direc filenames{i}],'/tractionNel');
                    f.Sn((i-1)*time_size+1:(i)*time_size,:) = temp;
                end


            case {'k','permeability'}
                % Find size to allocate
                temp = h5read([direc filenames{i}],'/permeability');
                xsize = size(temp,2);

                % Allocate
                f.k = zeros(time_size_full,xsize);

                % Read normal traction
                for i=1:length(filenames)
                    disp(['number of files loaded: ' num2str(i)])
                    temp = h5read([direc filenames{i}],'/permeability');
                    f.k((i-1)*time_size+1:(i)*time_size,:) = temp;
                end


            case {'stresstx','st'}
                % Find size to allocate
                temp = h5read([direc filenames{i}],'/tractionTel');
                xsize = size(temp,2);

                % Allocate
                f.St = zeros(time_size_full,xsize);

                % Read elastic shear traction
                for i=1:length(filenames)
                    disp(['number of files loaded: ' num2str(i)])
                    temp = h5read([direc filenames{i}],'/tractionTel');
                    f.St((i-1)*time_size+1:(i)*time_size,:) = temp;
                end

        end


    end
    return
    % By default read only V
else
    % Find size to allocate
    temp = h5read([direc filenames{i}],'/V');
    xsize = size(temp,2);

    % Allocate
    f.V = zeros(time_size_full,xsize);

    % Read V
    for i=1:length(filenames)
        disp(['number of files loaded: ' num2str(i)])
        temp = h5read([direc filenames{i}],'/V');

        f.V((i-1)*time_size+1:(i)*time_size,:) = temp;
    end
end
