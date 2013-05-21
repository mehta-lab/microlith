function img=computeimage(self,specimen,computeDevice)
% COMPUTEIMAGE: Computes images from the imaging system's transfer
% functions contained in the microlith object and the specimen properties.
%   IMG=COMPUTEIMAGE(MLOBJ, specimen,computeDevice) computes the image from the
%   transmission of the specimen/mask (in case of partially coherent
%   microscopy and lithography systems) or the density of fluorophores
%   (in case of incoherent microscopy) specified in specimen. The computed
%   image is returned and also saved within MLOBJ.
%
%   computeDevice = 'CPU', 'multiCPU, or 'GPU'. 
%   'multiCPU' option computes images across the focus in parallel if the
%   MATLAB parallel computing toolbox is installed and a matlabpool is
%   available.
% 
%   Current GPU implementation is only marginally faster than CPU and needs to be optimized.
% 
%   Written by Shalin Mehta, www.mshalin.com 
%   License: GPL v3 or later.

% This file is part of microlith package.
% 
%     microlith is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     microlith is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with microlith.  If not, see <http://www.gnu.org/licenses/>.

if(~isa(self,'microlith'))
    error('This function requires an object of the microlith class.');
end

L=length(self.x);

spec1d=isvector(specimen);

if(spec1d)
        % If the specimen is one dimensional vector, it is assumed constant
        % along the Y-axis.
        % For 1D specimen, the computation requires only 1D FFT, which
        % leads to substantial speed improvement for partially coherent
        % simulation.
        if(length(specimen)~=L)
            error('If the specimen is one dimensional, it has to be the same length as the axis over which imaging system is defined.');
        end
        
        if(iscolumn(specimen))
            specimen=specimen';
        end
        
        hsumOverSource=@sumOverSource1d; 
        % Function handle that calls 1D or 2D sum over source
        % implementation.
        %The sum over source implementation is rather different for 1D and 2D specimen, because oblique illumination along Y can affect the transmission of object spectrum along X.
        DCAlongY=floor(1+0.5*size(self.Po,1));
        ifftscale=L/numel(find(abs(self.Po(DCAlongY,:))>1E-20));

else
        if(size(specimen,1)~=L || size(specimen,1)~=L)
            error('If the specimen is two dimensional, it has has to be defined over the square grid.');
        end
        
        switch(computeDevice) % For 2D specimens, optimized functions that run on CPU and GPU are available.
            case {'CPU','multiCPU'}
                 hsumOverSource=@sumOverSource;
            case 'GPU'
                hsumOverSource=@sumOverSourceGPU;
            case 'CPUMex'
                hsumOverSource=@sumOverSource_mex;
        end
        
        % ifftscale for grid-independent calculation of image in coherent
        % and incoherent systems.
        ifftscale=L^2/numel(find(abs(self.Po)>1E-20));
        
     

end


        
if(spec1d)
        objspec=fftshift(fft(ifftshift(specimen)));
else
        objspec=fftshift(fft2(ifftshift(specimen)));
end

Tfun=self.Tfun;
img=zeros(size(specimen,1),size(specimen,2),size(Tfun,3),class(Tfun));
   % Radiometric factor that accounts for dependence of the image intensity
   % on imaging NA and size of the sensor pixel. The effect of the condenser
   % aperture relative to the imaging aperture is accounted for in the
   % sumOverSource function.
   % Our normalization is such that image of a 100nm^2 hole with an imaging
   % and illumination NA of 1 is unity.
   
   PixSize=self.x(2)-self.x(1);
   RadiometricFactor=(PixSize/0.1)^2*(self.params.NAo)^2;
   
%% Compute image.
switch(self.config)
    
    case {'Brightfield','Darkfield','PhaseContrast','DIC','DIC-Preza','PlasDIC'}
        if(strcmp(self.config,'DIC'))
            % Specimen spectrum is modulated according to shear and bias.
            shearangle=self.params.shearangle;
            halfshear=self.params.shear/2;
            halfbias=0.5*self.params.bias*pi/180;

            if(spec1d)
                    freqGridAlongX=self.m*cos(shearangle*pi/180);
                    objspec=1i*objspec.*sin(2*pi*freqGridAlongX*halfshear-halfbias); %TEST THIS.
            else
                    freqGridAlongShear=self.mm*cos(shearangle*pi/180)+...
                        self.nn*sin(shearangle*pi/180);
                    objspec=1i*objspec.*sin(2*pi*freqGridAlongShear*halfshear-halfbias);
            end
        end
        
        
   
                for idx=1:length(self.u)
                    img(:,:,idx)=RadiometricFactor*hsumOverSource(...
                        single(objspec),single(Tfun(:,:,idx)),single(self.Ic));
                end
%     
%             case 'GPU'
%                 objspecG=gpuArray(objspec);
%                 IcG=gpuArray(Ic);  
%                 for idx=1:length(self.u)
%                     TfunG=gpuArray(Tfun(:,:,idx));    
%                     imgG=hsumOverSource(objspecG,TfunG,IcG);
%                     img(:,:,idx)=gather(imgG);
%                 
                
   
        
    case {'Coherent','Fluorescence','Incoherent','Confocal'} %Tfun is CTF for coherent systems and OTF for incoherent systems.
        switch (computeDevice)
            case 'CPU'
                for idx=1:length(self.u)
                    img(:,:,idx)= RadiometricFactor*ifftscale*computeImageAmpInt(objspec,Tfun(:,:,idx)); 
                end
            case 'multiCPU'
                parfor idx=1:length(self.u)
                    img(:,:,idx)= RadiometricFactor*ifftscale*computeImageAmpInt(objspec,Tfun(:,:,idx)); 
                end
            case 'GPU'
                error('GPU optimization not yet implemented for coherent and incoherent systems. Please use multiCPU for speeding up computation.');
        end
        

end

%% Assign to the object.
switch(self.config)
    case {'Brightfield','Darkfield','PhaseContrast','DIC','DIC-Preza','PlasDIC','Coherent'}
        self.img=squeeze(img);
    case {'Fluorescence','Incoherent','Confocal'}
        self.img=squeeze(real(img));
        % Fluorescence image has to be real and therefore objspec has to be conjugate-symmetric.
        % But, I don't use 'symmetric option of ifft2 because of the
        % preceding ifftshift.
end
       


end

function img=computeImageAmpInt(objspec,Tfun) 
% Computes the inverse fourier transform of the product of object spectrum
% and imaging system's transfer function for coherent and incoherent
% imaging systems. Dimensionality of image is the same as dimensionality of
% the object. 
    if isvector(objspec)
            DCAlongY=floor(1+0.5*size(Tfun,1));
            img=fftshift(ifft(ifftshift(objspec.*Tfun(DCAlongY,:))));
    else
            img=fftshift(ifft2(ifftshift(objspec.*Tfun)));
    end
end
