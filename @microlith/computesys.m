function computesys(self,config,params)
% COMPUTESYS: Computes the transfer function of the imaging path and
% the intensity distribution of the illumination aperture for given
% microscopy/lithography system.    
%   COMPUTESYS(MLOBJ,config,params) calculates the transfer properties
%   using the config string and parameters specified in the params
%   structure.
%   
%   Valid string-values for config are:
%   'Brightfield','Darkfield','PhaseContrast','DIC','DIC-Preza','PlasDIC',
%   'Fluorescence','Incoherent','Confocal','Coherent', 'Custom'.
%   Following table summarizes the parameters tht can be present in the params structure 
%   and for which configurations they are used.
%     NAo           Numerical aperture of the imaging path (all configurations)
%     nImm          Refractive index of the objective immersion medium (all configurations)        
%     nEmbb         Refractive index of the specimen embedding medium (all configurations)
%     wavelength    Wavelength (in um) that dictates the resolution (all configurations)
%     NAc           Numerical aperture of the condenser (only for transmitted light methods).                           
%     annulus       For dark-field and phase-contrast systems the annulus
%                   parameter has to be a vector [InnerNA OuterNA].
%     phasering     For phase-contrast system phasering parameter has to be a
%                   vector [InnerNA OuterNA Absorption PhaseDelay]
%     shear         For DIC,DIC-Preza, and PlasDIC systems, 
%                   shear specified in units of wavelength/NAo.
%     bias          For DIC,DIC-Preza, and PlasDIC systems, 
%                   the phase-bias specified in degrees such that bias of
%                   90 leads to brightfield contrast.
%     shearangle    For DIC,DIC-Preza, and PlasDIC systems, 
%                   the direction of shear specified in degrees.
%     astigmatism   Astigmatism specified as coefficients of Zernike polynomials (n,m)=(2,-2) and
%     (2,2)
%     coma          Coma specified as coefficeients of Zernike polynomials
%     (n,m)= (3,1) and (3,-1).
%     spherical     Spherical aberration specified as coefficent for Zernike polynomial (n,m)=(4,0). 
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


%   The optical units used in our simulation are units in Born & Wolf  divided by 2pi. In
%   the units used by Born & Wolf, the first zero of airy disk occurs at
%   0.61*2pi=3.83 and the first on-axis minima occurs at normalized defocus of 4pi.
%   In the optical units used for simulation, the first zero of airy disk occurs at 0.61 and the
%   axial minima at 2.
%  
   % Establish normalized coordinates.
   %-----------------------------------------
   
   u=self.z*(params.NAo^2/params.wavelength)*(params.nEmbb/params.nImm^2);
   v=self.x*(params.NAo/params.wavelength);
   [self.vxx, self.vyy]=meshgrid(v);
   
    % Prepare the normalized spatial-frequency grid.
    % ---------------------------------------------------------
    L=length(v);
    vs=v(2)-v(1); 
    mc=1/(2*vs); % Spatial-frequency cutoff.
    self.m=linspace(-mc,mc,L); 
    self.m=self.m-min(abs(self.m)); % Ensure that DC position set to zero irrespective of numerical errors.
    [self.mm, self.nn]=meshgrid(self.m);
    
    Lu=length(u);
    if(Lu>1)
        us=u(2)-u(1);
        muc=1/(2*us);
        self.mu=linspace(-muc,muc,Lu);
    else 
        self.mu=NaN;
    end
    
    % Obtain azimuthal and radial co-ordinates in the pupil plane.
    %---------------------
    [mtheta, mr]=cart2pol(self.mm,self.nn);
    mrUnitCircle=mr.*(mr<=1); % Set all pixes outside of unit circle to zero for computation of Zernike polynomials.
   
   % Initialize pupils of the system
   %--------------------------------
    
    Po=single(mr<=1);
   % Po(Po<0)=1; %For pixels that are entirely within pupil set the transmittance to 1.
   % Po(Po>1)=0; %For pixels that are entirely outside pupil set the transmittance to 0. 
    
   % Model aberrations (astigmastism, coma, spherical) using Zernike
   % polynomials
   %---------------------------------
   if(isfield(params,'astigmatism'))
      Pastigm=params.astigmatism(1)*mrUnitCircle.^2.*cos(2*mtheta) +...
              params.astigmatism(2)*mrUnitCircle.^2.*sin(2*mtheta);  
   else
       Pastigm=zeros(size(Po));
   end
   
   if(isfield(params,'coma'))
       Pcoma=params.coma(1)*(3*mrUnitCircle.^3-2*mrUnitCircle).*cos(mtheta)+...
           params.coma(1)*(3*mrUnitCircle.^3-2*mrUnitCircle).*sin(mtheta);
   else
       Pcoma=zeros(size(Po));
   end
   
   if(isfield(params,'spherical'))
       Pspherical=params.spherical*(6*mrUnitCircle.^4-6*mrUnitCircle.^2+1);
   else
       Pspherical=zeros(size(Po));
   end
   
   PAberration=exp(1i*Pastigm).*exp(1i*Pcoma).*exp(1i*Pspherical);
   Po=Po.*single(PAberration);
  

     
    % Now set-up the pupils according to configuration
    %-----------------------
       
    Tfun=self.Tfun;
    switch(config)
            case {'Brightfield','DIC','DIC-Preza','PlasDIC'}
                S=params.NAc/params.NAo;
                Ic=single(mr<=S);  % Use of epsilon guards against numerical inaccuracy.             
                
    
                if(any(strcmp(config,{'DIC-Preza','PlasDIC'})))
                    
                    freqGridShear=self.mm*cos(params.shearangle*pi/180)+...
                        self.nn*sin(params.shearangle*pi/180);
                    halfshear=params.shear/2;
                    halfbias=0.5*params.bias*pi/180;
                    Po=1i*Po.*sin(2*pi*freqGridShear*halfshear-halfbias);                   
                
                    if(strcmp(config,'PlasDIC'))
                        PlasDICslit=abs(freqGridShear)<=(params.NAc/params.NAo);
                        % In PlasDIC, illumination is through a slit,
                        % instead of a circular aperture.
                        Ic=Ic.*double(PlasDICslit);
                    end
                    
                end
                

                for idx=1:length(u)
                    Tfun(:,:,idx)=Po.*exp(pi*1i*u(idx)*mr.^2);
                end
 

            case 'Darkfield'
                % Set the condenser annulus.
                if(~isfield(params,'annulus'))
                       error('For darkfield system, the params structure must have the ''annulus'' (a 1x2 array with innerNA and outerNA).');
                end
                
                Ic=single(mr>=params.annulus(1)/params.NAo ... Inner radius.
                       & mr<=params.annulus(2)/params.NAo ... Outer radius.
                       );
       
                parfor idx=1:length(u)
                    Tfun(:,:,idx)=Po.*exp(pi*1i*u(idx)*mr.^2);
                end

                
            case 'PhaseContrast'
                if(~ (isfield(params,'annulus') && isfield(params,'phasering')) )
                    error('For phase-contrast system, the params structure must have the ''annulus'' (1x2) and the ''phasering'' (1x4) vectors.');
                end
                % Set the condenser and objective annuli.
                % The pupil co-ordinate is normalized with respect to NA
                % of the objective.
                
                Ic=single(mr>=params.annulus(1)/params.NAo ... Inner radius.
                       & mr<=params.annulus(2)/params.NAo ... Outer radius.
                       );      
                
                phasering=params.phasering(3)* exp(1i*params.phasering(4)) *...
                single(mr>=params.phasering(1)/params.NAo & mr<=params.phasering(2)/params.NAo);
            
                %Direct multiplication by phasering sets transmission
                %outside phasering to zero. Therefore, create logical
                %mask to place the phasering in the objective pupil.
            
                phaseringmask=(phasering~=0);
                Po(phaseringmask)=phasering(phaseringmask);
                        
                parfor idx=1:length(u)
                    Tfun(:,:,idx)=Po.*exp(pi*1i*u(idx)*mr.^2);
                end
                
        case 'CustomPartiallyCoherent'%Replace Po and Ic by what was passed.
                Po=params.Po;
                Ic=params.Ic;
                parfor idx=1:length(u)
                    Tfun(:,:,idx)=Po.*exp(pi*1i*u(idx)*mr.^2);
                end
                                 
        case {'Fluorescence','Incoherent'}
            Ic=NaN;
            parfor idx=1:length(u)
                ctf=Po.*exp(pi*1i*u(idx)*mr.^2);
                Tfun(:,:,idx)=ctf2otf(ctf);
            end
            
        case {'CustomIncoherent'}
            Ic=NaN;
            Po=params.Po;
            parfor idx=1:length(u)
                ctf=Po.*exp(pi*1i*u(idx)*mr.^2);
                Tfun(:,:,idx)=ctf2otf(ctf);
            end
            
        case 'Confocal'
            Ic=NaN;
            parfor idx=1:length(u)
                ctf=Po.*exp(pi*1i*u(idx)*mr.^2);
                otfFL=ctf2otf(ctf);
                Tfun(:,:,idx)=ctf2otf(otfFL);
            end   
            
        case 'Coherent'
            Ic=double(self.mm==0 & self.nn==0);
            parfor idx=1:length(u)
                Tfun(:,:,idx)=Po.*exp(pi*1i*u(idx)*mr.^2);
            end
            
        case 'CustomCoherent'
            Ic=double(self.mm==0 & self.nn==0);
            Po=params.Po;
            parfor idx=1:length(u)
                Tfun(:,:,idx)=Po.*exp(pi*1i*u(idx)*mr.^2);
            end
            
        otherwise
            error(['The configuration: ' config ' is not implemented.']);
    end
    
    self.config=config;
    self.params=params;    
    self.Ic=Ic;
    self.Po=Po;
    self.Tfun=Tfun;
    self.u=u;
    self.v=v;
end
     
function otf=ctf2otf(ctf)
% OTF=CTF2OTF(CTF) computes the optical transfer function (OTF) from
% the coherent transfer function (CTF) of the imaging system.
% OTF is the auto-correlation of CTF, which we compute using the FFTs.
% The method generates the OTF on the same grid as CTF and hence the
% OTF has twice the cut-off frequency as that of CTF.
ctfSignificantPix=numel(find(abs(ctf)>eps(class(ctf))));
ifftscale=numel(ctf)/ctfSignificantPix;
%IFFT2 algorithm divides the input by numel(Po). The
%value at zero index of the output of IFFT2 is equal to the number of nonzero
%elements,i.e., numel(find(Po)). The above scale compensates
%for both and ensures that an image of a point produced by a clear
%circular aprture has a peak value of 1. This normalization allows us to
%compare images (apple to apple) computed over different grid sizes.

apsf=fftshift(ifft2(ifftshift(ctf)));
ipsf=ifftscale*abs(apsf).^2;

% The ifftscale should be applied to ipsf, because if applied to apsf; it
% is applied twice.

otf=fftshift(fft2(ifftshift(ipsf)));

end
