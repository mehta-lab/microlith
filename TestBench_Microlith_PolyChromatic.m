%% Evaluate accuracy of the Microlith package.
%
% 3D image of a point under fluorescent microscope with aberrations.
% Written by Shalin Mehta, www.mshalin.com, 
% License: GPL v3 or later.


%% Set-up simulation grid:

% All quantities are expressed in normalized optical coordinates 
% One can obtain the physical spatial coordinates like this:
% x (physical)= x (optical) * (wavelength/NA).
% In normalized coordinates, the jinc function's first zero occurs at 0.61
% and the circular pupil cuts-off at 1.

L=10; % Support over which we want to calculate the image. 

xs=0.05; % Sampling rate in the specimen plane.
% To avoid aliasing, specimen should be sampled at 0.2 wavelength/NA
% atleast.
% Sampling rate in the specimen plane determines extent in frequency
% domain.  xs=0.1 should suffice as the transfer function fits into the square with
% side 2 or (1+S)/sqrt(2) and sampling at 0.1 defines support of [-5 5].

v=-L:xs:L; % Transeverse extent of simulation.
u=-2:2*xs:2; % Axial extent of simulation.
[vx, vy]=meshgrid(v);

% Point specimen.
specimen=double(vx==0 & vy==0); 

% Some other specimens for which analytical image is easy to compute.
% % specimen=double(xx==0); %Slit
% % specimen=ones(size(xx)); %Transparent.


%% 3D PSF without chromatic aberration.

clear params;
params.NAo=1; 
params.wavelength=0.55; 
params.nImm=1.51; 
params.nEmbb=1.33;
RadioMetricFactor=params.NAo^2*(xs/0.1)^2;
% Factor that ensures radiometric consistency.

fluor=microlith(v,u);

computesys(fluor,'Fluorescence',params);
computeimage(fluor,specimen,'CPU');

%% 3D PSF with chromatic aberration.
params.wavelength=          [0.40    0.45    0.50    0.55    0.60    0.65    0.70];
params.CameraSensitivity=   [0.5     0.55    0.65    0.67    0.7     0.72    0.5];
params.SourceSpectrum=ones(size(params.wavelength)); % Source spectrum is irrelevant for widefield fluorescence imaging.
params.ChromaticShiftX=     [-0.2   -0.1     0        0      0       0.1     0.2];
params.ChromaticShiftY=zeros(size(params.wavelength));
params.ChromaticScaling=ones(size(params.wavelength));

fluorpoly=microlith(v,u);

computesys(fluorpoly,'Fluorescence',params);
computeimage(fluorpoly,specimen,'CPU');

%% Compare images.
figure(1); set(1,'Position',[100 100 1400 400],'color','white'); colormap hot;
for idx=10:30;
    clf;
    subplot(211);
    imagesc(v,v,fluor.img(:,:,idx),[0 RadioMetricFactor]); axis equal;
    title('Monochromatic'); xlim([-1 1]); ylim([-1 1]);
    
    subplot(212);
    imagesc(v,v,fluorpoly.img(:,:,idx),[0 RadioMetricFactor]); axis equal;
    title('With Chromatic Aberration'); xlim([-1 1]); ylim([-1 1]);
    pause(1);
end