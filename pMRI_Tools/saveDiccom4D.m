function [ mag, phase_1, phase_2, phase_3] = saveDiccom4D( directory_name, dirs, data)
%READDICCOM4D Summary of this function goes here
%   Detailed explanation goes here

% Create a new folder for the DICOM
mkdir('DICOM')
cd('DICOM')
mkdir(dirs.mag)
cd(dirs.mag)
save_dir.mag = pwd;
cd('..')

mkdir(dirs.phase_1)
cd(dirs.phase_1)
save_dir.phase_1 = pwd;
cd('..')

mkdir(dirs.phase_2)
cd(dirs.phase_2)
save_dir.phase_2 = pwd;
cd('..')

mkdir(dirs.phase_3)
cd(dirs.phase_3)
save_dir.phase_3 = pwd;
cd('../..')

cur_dir = pwd;
cd(directory_name)

%% Copy the Dicomdir file
copyfile('DICOMDIR',[save_dir.mag,'/../DICOMDIR'])

% Reshpae the data
sizes = size(data.mag);
data.mag = uint16(reshape(data.mag,[sizes(1), sizes(2), sizes(3)*sizes(4)]));
data.phase_1 = uint16(reshape(data.phase_1,[sizes(1), sizes(2), sizes(3)*sizes(4)]));
data.phase_2 = uint16(reshape(data.phase_2,[sizes(1), sizes(2), sizes(3)*sizes(4)]));
data.phase_3 = uint16(reshape(data.phase_3,[sizes(1), sizes(2), sizes(3)*sizes(4)]));

cd([directory_name,'/',dirs.mag])
saveImage(data.mag,save_dir.mag);

cd([directory_name,'/',dirs.phase_1])
saveImage(data.phase_1,save_dir.phase_1);

cd([directory_name,'/',dirs.phase_2])
saveImage(data.phase_2,save_dir.phase_2);

cd([directory_name,'/',dirs.phase_3])
saveImage(data.phase_3,save_dir.phase_3);

cd(cur_dir)


function saveImage(image,save_dir)
    dir_current = pwd;
    % Loop through the Diccom images 
    direct = dir(pwd);

    % Get som info from the Dicom
    num_images = length(direct) - 2;

    %% Import the data into a 4D array 

    % loop through all images in the directory and import them in the Instance
    % number order
    fprintf(sprintf('\nSaving data   %0.1f%%%%',0.0))
    for ind = 3:length(direct);
        
        % Get dicominfo
        info = dicominfo(direct(ind).name);  
        
%         image_no_order(:,:,info.InstanceNumber) = double(dicomread(direct(ind).name));
        cd(save_dir)
        dicomwrite(image(:,:,info.InstanceNumber), direct(ind).name, info, 'CreateMode', 'copy');
        cd(dir_current)
        % Display progress
        if ~mod(ind,20)
            if length(sprintf('%0.1f%%%%',100*ind/num_images))>5
                fprintf(sprintf('\b\b\b\b\b'))
                fprintf(sprintf('%0.1f%%%%',100*ind/num_images))
            else
                fprintf(sprintf('\b\b\b\b'))
                fprintf(sprintf('%0.1f%%%%',100*ind/num_images))
            end
        end
    end

end

end

