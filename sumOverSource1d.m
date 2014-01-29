function img=sumOverSource1d(objSpectrum,Po,srcInt,nshift,mshift,Ns,L)
% SUMOVERSOURCE1D Fast simulation of image of a 1D specimen under partially
% coherent system.
%   IMG=SUMOVERSOURCE1D(objSpectrum,Po,srcInt,nshift,mshift,Ns,L) computes partially
%   coherent image of 1D specimen, objSpectrum is a 1D vector instead of 2D matrix
% 
% 
%   See also: sumOverSource.
%
%   Written by Shalin Mehta, www.mshalin.com 
%   License: GPL v3 or later.

    
% Algorithm for computing image of 1D specimen. 
% Depending on position of source point in Y, identify the line of Po that
% affects the image.
% Depending on position of source point in X, shift object spectrum along X.
% Multiply and perform 1D inverse FFT to compute image.

    % I assume that all the matrices are defined over the same frequency grid.
    % This function is not immune to errors. Error checks must be
    % implemented by the parent function. The error checks are ommited
    % so that the function can be run on GPU without having to transfer any
    % data to CPU.
    objspecshifted=objSpectrum;
    objspecshifted(:)=0;
    img=objSpectrum;
    img(:)=0;
    
    % The spectral support (or spatial cutoff) of simulation must be able
    % to accommodate the largest spatial frequency that is captured.
    % Therefore, there is no need to pad either objspec or Po to find the
    % product of shifted objspec and Po.
        
    % For each source point, calculate what indices from the unshifted
    % spectrum map to what indices of shifted spectrum.
    %---------------------------------------
    
    % Where do shifted and unshifted spectra start given the sign of mshift.
    %---------------------------------------------------------------
    negmshift=single(mshift<0); posmshift=single(mshift>=0);
    
    mSpectrumStart=negmshift.*(-mshift+1) +...
        posmshift.*1;
    mSpectrumEnd=negmshift.*(L)+ ...
        posmshift.*(L-mshift);
    
    mShiftSpectrumStart=negmshift.*1 +...
        posmshift.*(mshift+1);
    mShiftSpectrumEnd=negmshift.*(L+mshift) +...
        posmshift.*L;
        
%     % Above is the vectorized form of the following logic for each mshift.
%     % It is vectorized to run efficiently on GPU.
%          if(mshift<0)
%              mSpectrumidx= (-mshift+1):L;
%              mShiftedSpectrumidx= 1:(L+mshift);
% 
%          else
%              mSpectrumidx=1:(Len-mshift);
%              mShiftedSpectrumidx=(1+mshift):Len;
%          end
    
    loopvar=ones(1,'like',Ns);

    
    while loopvar<=Ns
        % A parallel computing demo at http://www.mathworks.com/products/demos/parallel-computing/paralleldemo_gpu_mandelbrot/
        % suggests that while loop is more efficient on GPU then for
        % loop.

         nidxPo=nshift(loopvar)+floor(0.5*L); %nshift of zero corresponds to the center of Po.
         
         midxshifted=colon(mShiftSpectrumStart(loopvar),mShiftSpectrumEnd(loopvar));         
         midxunshifted=colon(mSpectrumStart(loopvar),mSpectrumEnd(loopvar));
         
         objspecshifted(midxshifted)=objSpectrum(midxunshifted);
         
         imgspectrum=objspecshifted.*Po(nidxPo,:);
           

           %IFFT2 algorithm divides the input by numel(Po). The
           %value at zero index of the output of IFFT2 is equal to the number of nonzero
           %elements,i.e., numel(find(Po)). The above scale compensates
           %for both and ensures that an image of a point produced by a clear
           %circular aprture has a peak value of 1. This normalization allows us to
           %compare images (apple to apple) computed over different grid sizes.


           imgpoint=fftshift(ifft(ifftshift(imgspectrum ))); 
           imgpoint=srcInt(loopvar)*abs(imgpoint).^2;
           img=img+imgpoint;
           
           % As we sum over the source-points, the peak value in img depends on the sampling of
           % source. In the end, we normalize by total number of transmissive pixels in the
           % objective pupil to ensure that the image intensity depends
           % only on the shape of the condenser pupil and not on the
           % sampling interval. This is important when we want to compare
           % the images computed on different grid sizes.
           
%           debug=2; 
%               figure(debug);
%                %xc=1/(2*dm); xnorm=linspace(-xc,xc,length(mnorm));
%                subplot(221); imagesc(abs(objspecshifted)); axis equal; colorbar;
%                title(['Shifted objSpectrum: (\mm \nn)=' num2str([mshift(loopvar) nshift(loopvar)])]);
%                subplot(222); imagesc(log(1+abs(imgspectrum))); axis equal; colorbar;
%                title('Filtered objSpectrum');
%                subplot(223); imagesc(abs(imgpoint).^2); axis equal;  colorbar;
%                %xlim([-Len Len]); ylim([-Len Len]);
%                title('Image intensity due to this condenser point');
%                subplot(224); imagesc(img); axis equal; colorbar;
%                title('Sum image upto now');
%                %xlim([-Len Len]); ylim([-Len Len]);
%                pause(0.01);
          
          loopvar=loopvar+1;
        
    end
    
% Normalize the image such that the values do not depend on the fineness of
% the source grid.

% Following is the normalization according to Martin's book. We do not use it, because the
% following leads to divide by zero for dark-field system.
%         normfactor=abs(Po).^2.*abs(Ic); 
%         normfactor=sum(normfactor(:));
%         if(normfactor) %If normafactor == 0, we are imaging with dark-field system.
%             img=img/sum(normfactor(:));
%         end

end
