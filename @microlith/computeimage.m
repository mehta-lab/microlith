function img=computeimage(self,specimen,computeDevice)
% COMPUTEIMAGE: Computes images from the imaging system's transfer
% functions contained in the microlith object and the specimen properties.
%   IMG=COMPUTEIMAGE(MLOBJ, specimen,computeDevice) computes the image from the
%   transmission of the specimen/mask (in case of partially coherent
%   microscopy and lithography systems) or the density of fluorophores
%   (in case of incoherent microscopy) specified in specimen. The computed
%   image is returned and also saved within MLOBJ.
%
%   computeDevice = 'CPU', 'multiCPU', or 'GPU'. 
%   'multiCPU' option computes images across the focus in parallel if the
%   MATLAB parallel computing toolbox is installed and a matlabpool is
%   available.
% 
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

% Size and center of grid.
L=length(self.x);
DCAlongY=floor(1+0.5*size(self.Po,1));
DCAlongX=floor(1+0.5*size(self.Po,2));


% Is specimen 1D?
spec1d=isvector(specimen);

% Obtain the specimen spectrum. The pupils are stored in microlith object.
% Choose the right image computation algorithm.
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
        
        objspec=fftshift(fft(ifftshift(specimen)));
        
        hsumOverSource=@sumOverSource1d; 
        % Function handle that calls 1D or 2D sum over source
        % implementation.
        %The sum over source implementation is rather different for 1D and 2D specimen, because oblique illumination along Y can affect the transmission of object spectrum along X.
        

        
else
        if(size(specimen,1)~=L || size(specimen,1)~=L)
            error('If the specimen is two dimensional, it has has to be defined over the square grid.');
        end
        
        objspec=fftshift(fft2(ifftshift(specimen)));
         hsumOverSource=@sumOverSourcePar;

end

% Read the pupils and allocate memory for the image.
%---------------------------------------------------
% We use single precision representation, because that runs the fastest on
% GPU. The round-off errors due to single-precision are not too expensive.

Tfun=single(self.Tfun);
img=zeros(size(specimen,1),size(specimen,2),size(Tfun,3),'single');

% Configuration specific adjustments.
%----------------------------------
% Specimen spectrum is modulated according to shear and bias in DIC.
if(strcmp(self.config,'DIC'))
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

   
%% Compute image.
switch(self.config)
    
    case {'Brightfield','Darkfield','PhaseContrast','DIC','DIC-Preza','PlasDIC','CustomPartiallyCoherent'}

            %Find non-zero pixels on condenser and the amount of spectral shift
            %induced by them.
            %------------------------------------------
            Ic=single(self.Ic);
            [nshift,mshift, srcInt]=find(Ic);
            nshift=nshift-DCAlongY;
            mshift=mshift-DCAlongX;
            
            % Use single precision floating point to improve speed/memory.
            specimen=single(specimen);
            objspec=single(objspec);
            nshift=single(nshift);
            mshift=single(mshift);
            Ns=single(numel(srcInt));
            switch(computeDevice)
                case 'CPU' % Serial computation. FFT is multi-threaded.
                    parfor idx=1:length(self.u)
                        %img(:,:,idx)=sumOverSource(objspec,Tfun(:,:,idx),Ic);
                        img(:,:,idx)=hsumOverSource(objspec,Tfun(:,:,idx),srcInt,nshift,mshift,Ns,L);
                    end          
                case 'GPU' % Parallel computation with GPU.
                    % Transfer the vectors and matrices to GPU.
                    objspecG=gpuArray(objspec);
                    nshiftG=gpuArray(nshift); 
                    mshiftG=gpuArray(mshift);
                    srcIntG=gpuArray(srcInt);
                    NsG=gpuArray(Ns);
                    Lg=gpuArray(L);
                    for idx=1:length(self.u)
                        PoG=gpuArray(Tfun(:,:,idx));
                        %imgG=hsumOverSource(objspecG,PoG,srcIntG,nshiftG,mshiftG,NsG,Lg);
                        imgG=sumOverSourcePar(objspecG,PoG,srcIntG,nshiftG,mshiftG,NsG,Lg);
                        img(:,:,idx)=gather(imgG);
                        clear('PoG'); % Clear to make space on GPU.
                    end
                    clear('objspecG','nshiftG','nshiftG','mshiftG','srcIntG');
%                 case 'CPUvectorized'
%                     for idx=1:length(self.u)
%                         img(:,:,idx)=sumOverSourceVectorized(specimen,Tfun(:,:,idx),srcInt,nshift,mshift,Ns);
%                         %img(:,:,idx)=sumOverSourceUni(objspec,Tfun(:,:,idx),srcInt,nshift,mshift,Ns);
%                     end                       

            end
  
    case {'Coherent','Fluorescence','Incoherent','Confocal','CustomCoherent','CustomIncoherent'} %Tfun is CTF for coherent systems and OTF for incoherent systems.
        switch (computeDevice)
            case 'CPU'
                for idx=1:length(self.u)
                    img(:,:,idx)= computeImageAmpInt(objspec,Tfun(:,:,idx)); 
                end
            case 'multiCPU'
                parfor idx=1:length(self.u)
                    img(:,:,idx)= computeImageAmpInt(objspec,Tfun(:,:,idx)); 
                end
            case 'GPU'
                error('GPU optimization not yet implemented for coherent and incoherent systems. Please use multiCPU for speeding up computation.');
        end
end

% Apply the radiometric factor and assign to the object.
%-------------------------------------------------

% Radiometric factor that accounts for dependence of the image intensity
% on imaging NA and size of the sensor pixel. The effect of the condenser
% aperture relative to the imaging aperture is accounted for in the
% sumOverSource function.
% Our normalization is such that image of a 100nm^2 hole with an imaging
% and illumination NA of 1 is unity. For any other illumination NA, the
% intensity is proportional to ratio of the area of the condenser pupil to
% area of the objective pupil.

if(spec1d)
        % ifftscale for grid-independent calculation of image in coherent
        % and incoherent systems.
        Np=numel(find(abs(self.Po(DCAlongY,:))>1E-12));
        ifftscale=L/Np;   

else
        % ifftscale for grid-independent calculation of image in coherent
        % and incoherent systems.
        Np=numel(find(abs(self.Po)>1E-12));
        ifftscale=L^2/Np;       
           %IFFT2 algorithm divides the input by numel(Po). The
           %value at zero index of the output of IFFT2 is equal to the number of nonzero
           %elements,i.e., numel(find(Po)). The above scale compensates
           %for both and ensures that an image of a point produced by a clear
           %circular aprture has a peak value of 1. This normalization allows us to
           %compare images (apple to apple) computed over different grid sizes.
           
           % Note that the numel(find(Po)) should be substituted by
           % numel(find(Po*ObjectSpectrum)) when the object spectrum is
           % finite and when the objectspectrum is shifted during partially coherent computation.
           % However, for realistic specimens, the object spectrum almost
           % always exceeds twice the support of Po and therefore the
           % product of the pupil and shifted spectrum has the same support
           % as the pupil.


end


PixSize=self.x(2)-self.x(1);


switch(self.config)
    case 'Coherent'
        RadiometricFactor=ifftscale*(PixSize/0.1)^2*(self.params.NAo)^2; 
        img=squeeze(RadiometricFactor*img);  

    case {'Brightfield','Darkfield','PhaseContrast','DIC','DIC-Preza','PlasDIC'}
        Sfactor=1/Np; % sum over soure 'brightens' the image by number of source points. 
        % Normalization by the number of transparent points in imaging
        % pupil ensures that the intensity is proportional to Ns/Np.        
        RadiometricFactor=ifftscale^2*Sfactor*(PixSize/0.1)^2*(self.params.NAo)^2; 
        % For partially coherent imaging, ifft is applied on amplitude
        % image and then the result is squared. The radiometric factor
        % needs to take into account this squaring.
        img=squeeze(RadiometricFactor*img);
        
    case {'Fluorescence','Incoherent','Confocal'}
        RadiometricFactor=ifftscale*(PixSize/0.1)^2*(self.params.NAo)^2; 
        img=squeeze(real(RadiometricFactor*img));        
        % Fluorescence image has to be real and therefore objspec has to be conjugate-symmetric.
        % But, I don't use 'symmetric option of ifft2 because of the
        % preceding ifftshift.
end
    
%%%% If polychromatic illumination.
%%---------------------------------

if(~isscalar(self.params.wavelength)) 
    SpectralWeights=gray2norm(self.params.SourceSpectrum).*gray2norm(self.params.CameraSensitivity);

    CenterWavelengthIdx=find(self.params.wavelength == self.params.CentralWavelength);
    
    SpectralWeights=SpectralWeights/sum(SpectralWeights); % Ensure area unde spectral weights is unity.
    xshifts=self.params.ChromaticShiftX-self.params.ChromaticShiftX(CenterWavelengthIdx);
    yshifts=self.params.ChromaticShiftY-self.params.ChromaticShiftY(CenterWavelengthIdx);
    DiffractionScaling=self.params.wavelength/self.params.CentralWavelength;
    ChromaticScaling=self.params.ChromaticScaling/self.params.ChromaticScaling(CenterWavelengthIdx);
    Scaling=DiffractionScaling.*ChromaticScaling;
    for idw=1:numel(self.params.wavelength)
        affinemat=ShiftScaleRotToaffine(xshifts(idw),yshifts(idw),Scaling(idw),Scaling(idw),0);
        imgPresentWavelength=SpectralWeights(idw)*imtransformAffineMat(double(img),affinemat,'bilinear','coordinates','centered','PixSize',PixSize);
        self.img=self.img+imgPresentWavelength;
    end
else
    self.img=img;
    
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
