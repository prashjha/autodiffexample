% A script which calculates a relative flip angle map of the 13C clamshell coil
% given specific scan prameters from clinical HP MRI patient scans. 
% 
% Collin Harlan
% 04/06/2023

clear
close
clc

%% Load in eight "layers" of 2D matrices (eight measurment points in y direction) hand measured from the clamshell coil

for ii=1:8
    eval(sprintf('load(''layer_%d_modified.mat'');', ii))
end

%% Convert measured values from dB to B1

layer_1_modified = (10.^(layer_1_modified/(20)));
layer_2_modified = (10.^(layer_2_modified/(20)));
layer_3_modified = (10.^(layer_3_modified/(20)));
layer_4_modified = (10.^(layer_4_modified/(20)));
layer_5_modified = (10.^(layer_5_modified/(20)));
layer_6_modified = (10.^(layer_6_modified/(20)));
layer_7_modified = (10.^(layer_7_modified/(20)));
layer_8_modified = (10.^(layer_8_modified/(20))); 

%% Concatenate the eight layers to create a single 3D matrix 

dat = cat(3, layer_1_modified, layer_2_modified, layer_3_modified, layer_4_modified, layer_5_modified, layer_6_modified, layer_7_modified, layer_8_modified);

%% Create a map of the clamshell volume over which S21 measurments were taken in x,y, and z directions 

x = linspace(-14*2.54/2, 14*2.54/2, 14); %14 points measured over 14 inches, convert to cm with 2.54 conversion factor
y = ([1.5875, 2.9625, 4.15, 5.3375, 6.525, 7.7125, 8.9, 10.0875] - 5.25) * 2.54; % eight measurments were taken in the y direction, convert to cm with 2.54 conversion factor
z = linspace(-17*2.54/2, 17*2.54/2, 17); %17 points measured over 17 inches, convert to cm with 2.54 conversion factor
[x_grid, y_grid, z_grid] = meshgrid(x, y, z);

%% Normalize B1 data to interpolated point value at the center (-1.5, 0.22, -1.5) of the 3D matrix created above 

point_val = interp3(x_grid, y_grid, z_grid, permute(dat,[3,2,1]), -1.5, 0.22, -1.5);
dat_norm = 100 * dat ./ point_val;
dat_norm = permute(dat_norm,[3,2,1]);

%% Generate a map showing flip angle changes across the FOV 

% Givens: we have a 16x16 matrix, 24cm FOV, and 8 slices that are 1.5cm
% thick. Assume that we are operating at isocenter.

% Create 3D volumetric matrix map using given scan parameters
FOV = 24; %24 cm FOV for x and y direction
matrix_size = 16; %16x16 matrix
pixel_size = FOV/matrix_size; %cm^2
slices = 8;
slice_thickness = 1.5; %cm
axial_distance = slices*slice_thickness; %12 cm for z direction

x_new = linspace(-FOV/2, FOV/2, matrix_size); %sagittal, yz plane
y_new = linspace(-FOV/2, FOV/2, matrix_size); %coronal, xz plane
z_new = linspace(-axial_distance/2, axial_distance/2, slices); %axial, xy plane

[x_grid_new, y_grid_new, z_grid_new] = meshgrid(x_new, y_new, z_new);

%% interp original points to new 3D volume created using scan parameters

% This is a relative flip angle map!
Flip_Angle_Map = interp3(x_grid, y_grid, z_grid, dat_norm, x_grid_new, y_grid_new, z_grid_new)

pyruvateFA = 20.*.01*Flip_Angle_Map;
lactateFA = 30.*.01*Flip_Angle_Map;
niftiwrite(pyruvateFA,'famappyr.nii')
niftiwrite(lactateFA ,'famaplac.nii')
min(pyruvateFA(:))
max(pyruvateFA(:))

% NOTE: The measurements of the clamshell are as follows:
% 
% X = 35.56 cm
% Y = 22 cm
% Z = 43.18 cm
% 
% Becuase this is an axial scan (z direction into the bore), the 3D volume given the scan parameters is 
% 
% X = 24 cm (FOV)
% Y = 24 cm (FOV)
% Z = 12 cm (slices*thickness)
% 
% So when you do the mapping, you get NaN’s for the first two rows because 24 cm > 22 cm.

% This should not affect the analysis because we don’t typically have patient 
% anatomy or tumors of interest at the edge of the FOV.
















