function img=sumOverSourceGPU(objSpectrum,Po,Ic)
% sumOverSourceGPU GPU implementation of the sumOverSource algorithm. Needs
% to be re-written, any inputs from users welcome. 
% Currently, it is only about a factor of 4 faster than CPU implementation. 
%
%   See also: sumOverSource.
%
%   Written by Shalin Mehta, www.mshalin.com 
%   License: GPL v3 or later.
    
    % I assume that all the matrices are defined over the same frequency grid.
    % This function is not immune to errors. Error checks must be
    % implemented by the parent function. The error checks are ommited
    % so that the function can be run on GPU without having to transfer any
    % data to CPU.
    
    %% Transfer threee matrices to GPU,
    objspecG=gpuArray(objSpectrum);
    IcG=gpuArray(Ic);  
    PoG=gpuArray(Po);

    %% Do all the computation on GPU.
    %--------------------------------
    objspecshifted=PoG;
    objspecshifted(:)=0;
    imgG=objspecG;
    imgG(:)=0;
    
    % The spectral support (or spatial cutoff) of simulation must be able
    % to accommodate the largest spatial frequency that is captured.
    % Therefore, there is no need to pad either objspecG or PoG to find the
    % product of shifted objspecG and PoG.
    
    %Find non-zero pixels on condenser and the amount of spectral shift
    %induced by them.
    %------------------------------------------
    [nshift,mshift, SrcInt]=find(IcG);
    
    FreqGridLen=gpuArray(size(IcG,1));    
    freqOrigin=floor(0.5*FreqGridLen+1);
    nshift=nshift-freqOrigin;
    mshift=mshift-freqOrigin;
    
    % For each source point, calculate what indices from the unshifted
    % spectrum map to what indices of shifted spectrum.
    %---------------------------------------
    
    % Where do shifted and unshifted spectra start given the sign of mshift.
    %---------------------------------------------------------------
    negmshift=single(mshift<0); posmshift=single(mshift>=0);
    
    mSpectrumStart=negmshift.*(-mshift+1) +...
        posmshift.*1;
    mSpectrumEnd=negmshift.*(FreqGridLen)+ ...
        posmshift.*(FreqGridLen-mshift);
    
    mShiftSpectrumStart=negmshift.*1 +...
        posmshift.*(mshift+1);
    mShiftSpectrumEnd=negmshift.*(FreqGridLen+mshift) +...
        posmshift.*FreqGridLen;
        
%     % Above is the vectorized form of the following logic for each mshift.
%     % It is vectorized to run efficiently on GPU.
%          if(mshift<0)
%              mSpectrumidx= (-mshift+1):FreqGridLen;
%              mShiftedSpectrumidx= 1:(FreqGridLen+mshift);
% 
%          else
%              mSpectrumidx=1:(Len-mshift);
%              mShiftedSpectrumidx=(1+mshift):Len;
%          end
    
      % Likewise for nshift
      %--------------------------
      negnshift=single(nshift<0); posnshift=single(nshift>=0);
      nSpectrumStart=negnshift.*(-nshift+1) +...
        posnshift.*1;
      nSpectrumEnd=negnshift.*(FreqGridLen)+ ...
        posnshift.*(FreqGridLen-nshift);
    
      nShiftSpectrumStart=negnshift.*1 +...
        posnshift.*(nshift+1);
      nShiftSpectrumEnd=negnshift.*(FreqGridLen+nshift) +...
        posnshift.*FreqGridLen;
    
        loopvar=gpuArray(1);
        Ns=gpuArray(numel(nshift));


    %waith=waitbar(loopvar/Ns,['Computing image intensity by summing over the source:' num2str(100*loopvar/Ns,2) '%']);
    while loopvar<=Ns
        % A parallel computing demo at http://www.mathworks.com/products/demos/parallel-computing/paralleldemo_gpu_mandelbrot/
        % suggests that while loop is more efficient on GPU then for
        % loop.

        %waitbar(loopvar/Ns,waith,['Computing image intensity by summing over the source:' num2str(100*loopvar/Ns,2) '%']);
         nidxshifted=gpuArray.colon(nShiftSpectrumStart(loopvar),nShiftSpectrumEnd(loopvar));
         midxshifted=gpuArray.colon(mShiftSpectrumStart(loopvar),mShiftSpectrumEnd(loopvar));
         
         nidxunshifted=gpuArray.colon(nSpectrumStart(loopvar),nSpectrumEnd(loopvar));
         midxunshifted=gpuArray.colon(mSpectrumStart(loopvar),mSpectrumEnd(loopvar));
         
         objspecshifted(nidxshifted,midxshifted)=objspecG(nidxunshifted,midxunshifted);
         
         imgspectrum=objspecshifted.*PoG;
           
         imgspectrumsignificant=abs(imgspectrum)>1E-20;
         ifftscale=FreqGridLen^2/sum(imgspectrumsignificant(:));
          
           %IFFT2 algorithm divides the input by numel(PoG). The
           %value at zero index of the output of IFFT2 is equal to the number of nonzero
           %elements,i.e., numel(find(PoG)). The above scale compensates
           %for both and ensures that an image of a point produced by a clear
           %circular aprture has a peak value of 1. This normalization allows us to
           %compare images (apple to apple) computed over different grid sizes.


           imgpoint=fftshift(ifft2(ifftshift( ifftscale*imgspectrum ))); 
           imgpoint=SrcInt(loopvar)*abs(imgpoint).^2;
           imgG=imgG+imgpoint;
           
           % As we sum over the source-points, the peak value in imgG depends on the sampling of
           % source. In the end, we normalize by total number of transmissive pixels in the
           % objective pupil to ensure that the image intensity depends
           % only on the shape of the condenser pupil and not on the
           % sampling interval. This is important when we want to compare
           % the images computed on different grid sizes.
           
%           debug=2; 
%               figure(debug);
%                %xc=1/(2*dm); xnorm=linspace(-xc,xc,length(mnorm));
%                subplot(221); imagesc(abs(objspecshifted)); axis equal; colorbar;
%                title(['Shifted objspecG: (\mm \nn)=' num2str([mshift(loopvar) nshift(loopvar)])]);
%                subplot(222); imagesc(log(1+abs(imgspectrum))); axis equal; colorbar;
%                title('Filtered objspecG');
%                subplot(223); imagesc(abs(imgpoint).^2); axis equal;  colorbar;
%                %xlim([-Len Len]); ylim([-Len Len]);
%                title('Image intensity due to this condenser point');
%                subplot(224); imagesc(imgG); axis equal; colorbar;
%                title('Sum image upto now');
%                %xlim([-Len Len]); ylim([-Len Len]);
%                pause(0.01);
          
          loopvar=loopvar+1;
        
    end

  %  delete(waith);
% Normalize the image such that the values do not depend on the fineness of
% the source grid.

    imgG=imgG/Ns;

    %% Get the result from GPU and return to the calling function.
    img=gather(imgG);
% Following is the normalization according to Martin's book. We do not use it, because the
% following leads to divide by zero for dark-field system.
%         normfactor=abs(PoG).^2.*abs(IcG); 
%         normfactor=sum(normfactor(:));
%         if(normfactor) %If normafactor == 0, we are imaging with dark-field system.
%             imgG=imgG/sum(normfactor(:));
%         end

end
