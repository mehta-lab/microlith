function mutint=computecoherence(self)
% COMPUTECOHERENCE: Computes the mutual intensity at the focal region of
% the microscope.
%   mutint=COMPUTECOHERENCE(MLOBJ) computes the mutual intensity (complex coherence)
% from the intensity distribution of the source, by employing the idea of
% gernealized source developed by C. W. McCutchen.
%
%  Reference: C. McCUTCHEN, "Generalized Source and the van Cittert?Zernike
%  Theorem: A Study of the Spatial Coherence Required for Interferometry,"
%  J. Opt. Soc. Am.  56, 727-732 (1966).
%
%  The computed mutual intensity is returned as complex-valued matrix.
%
%   Written by Shalin Mehta, www.mshalin.com 
%   License: GPL v3 or  later.

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


switch(self.config)
    
    case {'Brightfield','Darkfield','PhaseContrast','DIC','KohlerDIC','PlasDIC','Coherent'}
        
        Ic=self.Ic;
        
        if(strcmp(self.config,'DIC'))
            % For DIC mutual coherence is between two polarized beams. It can be
            % visualized by modulating condenser pupil with sinusoid
            % depending on the shear and bias.
        end
        
        L=length(self.x);
        mutint=zeros(L,L,length(self.u));

        ifftscale=numel(self.Ic)/numel(find(abs(self.Ic)>1E-20));
        % Obtain azimuthal and radial co-ordinates in the pupil plane.
        [mtheta, mr]=cart2pol(self.mm,self.nn);%#ok<ASGLU>

        u=self.u;
        parfor idx=1:length(u)
            GenSource=Ic.*exp(pi*1i*u(idx)*mr.^2);
            mutint(:,:,idx)=fftshift(ifft2(ifftshift(ifftscale*GenSource)));
        end

        
    case {'Fluorescence','Incoherent','Confocal'}
        error('Calculation of mutual cohrence is not applicable to incoherent imaging systems.');

end


end
