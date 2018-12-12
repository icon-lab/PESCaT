%% PESCaT Demo
% Reconstruction by Projetion over Epigraph Sets and Calibration over
% Tensors (PESCaT) Demo
%
% This is a demo that shows an example implementation
% of PESCaT algorithm over bSSFP data with multiple phase-cycled
% acquisitions and coils. It is based on Shahdloo et. al, "Projection onto
% Epigraph Sets for Rapid Self-Tuning Compressed Sensing MRI". PESCaT is a
% technique that employs projection onto epigraph sets of L1 and TV norm
% functions to self-tune the regularization parameters. A tensor
% interpolation kernel is estimated over the aggregated data across coils
% and phase cycles and used to fill missing k-space points.

%%

clearvars;
close all;

addpath(genpath('ESPIRiT'))
%% Parameters
% The provided dataset in this tutorial have four coils and eight datasets.
% In multi-acquisition bSSFP imaging, each phase-cycled acquisition is
% accelerated by a factor equal to the total number of acquisitions to keep
% the acquisition time intact. Here we use N=4 phase-cycles and R=4
% acceleration factor. In the end, images are p-norm combined with p=2
% across coils and p=4 across acquisitions.
% To compare the PESCaT reconstruction with convensional reconstructions
% using the SURE criterion, and alternative reconstruction is also
% performed.

N = 4; % number of phase-cycled acquisitions
R = 4; % acceleration factor
p_acq = 4; % p value for p-norm combination over phase-cycles
p_coils = 2; % p value for p-norm combination over coils

reconSURE = 1; % do you want SURE reconstruction too?

%% Necessary libraries and folders
% SPIRiT V0.3 is required for the wavelet operations, check for it.
if exist('crop.m', 'file')==0
	addpath(genpath('ESPIRiT'));
	if exist('crop.m', 'file')==0
		warning('ReCat utilizes SPIRiT V0.3 library. Please download from http://people.eecs.berkeley.edu/~mlustig/Software.html and copy the "ESPIRiT" directory to the same directory as this demo file. ');
		return;
	end
end
%%
% Adding the util folder that includes some utility functions
if exist('normalize.m', 'file')==0
	addpath('util');
end

%% Loading tutorial data
% load sample bSSFP data
load('data/invivo_4coil.mat');
raw_data = double(raw_data); % LSQR implementation requires double type

%% Loading undersampling Mask
% load pre-generated mask
load(['masks/mask_' num2str(R) 'x.mat']);

%% Reference Image
% Generate the fully-sampled reference via p-norm combination over N=4 of
% the fully-sampled acquisitions

images = ifft2c(raw_data(:,:,1:2:8,:));
originalImage = sos(sos(images,4,p_coils),3,p_acq);

%% Prepare Data
% Take N-many acquisitions
imageFFT = reshape(raw_data(:,:,1:8/R:8,:),[size(images,1),size(images,2),1,N,size(images,4)]);
%%
% In reality, the sampling masks are different across acquisitions but the
% same across slices and coils. Undersampled data is assessed by
% multiplying the k-space data with the undersampling mask.
mask = repmat(mask(:,:,1:8/R:8),[1,1,1,size(imageFFT,3),size(imageFFT,5)]);
mask = permute(mask,[1,2,4,3,5]);
sampling.mask = mask;
maskedData = imageFFT.*mask;

%% Create PESCaT instance
% A PESCaT object is instantiated with the default properties. These can be
% modified by passing extra arguments to the object instantiator.
pobj = PESCaT(maskedData,sampling);

%% Perfrorm reconstruction
% The reconPESCaT method performs the reconstruction on the data passed to
% the PESCaT object.
pobj.reconPESCaT();

%% Outputting the results
% Normalize both reference image and reconstructed image to 0-1.
originalImage = normalize(originalImage);
result = normalize(pobj.recon);
%%
% Display the reconstruction and the reference images.
figure; imshow(originalImage);
title('Fully Sampled Image');
figure; imshow(result);
title('PESCaT Reconstructed Image');
%%
% Print the PSNR using the fully-sampled reference.
fprintf('PESCaT PSNR: %.2f\n', psnr(result, originalImage));
%%
% Print the elapsed time.
fprintf('PESCaT elapsed time: %.2f\n', sum(pobj.optimParams.elapsed));

%% Alternative reconstruction
% This conditional statement performs the reconstruction using the
% convensional method and outputs the results.

if reconSURE
    %%
    % Create PESCaT instance with sparsity and TV projections set to `SURE' and 'STD'
    sobj = PESCaT(maskedData,sampling, 'sparsityType', 'SURE', 'TVType', 'STD');
    %%
    % Perfrorm reconstruction
    sobj.reconPESCaT();
    resultSURE = normalize(sobj.recon);
    %%
    % Output the results
    figure; imshow(originalImage);
    title('Fully Sampled Image');
    figure; imshow(resultSURE);
    title('SURE Reconstructed Image');
    
    fprintf('SURE PSNR: %.2f\n', psnr(resultSURE, originalImage));
    fprintf('SURE elapsed time: %.2f\n', sum(sobj.optimParams.elapsed));
end
