function img=sumOverSourcePar(objSpectrum,Po,srcInt,nshift,mshift,Ns,L)
% sumOverSourceUni implements easily parallelized version of sum over
% source algorithm. This modified algorithm can be run on GPU and
% multi-core CPU without incurring excessive memory transfer overheads.
% If the input arguments are gpuArrays, the algorithm runs on GPU;
% otherwise it runs on CPU.
%   See also: sumOverSource.
%
%   Written by Shalin Mehta, www.mshalin.com 
%   License: GPL v3 or later.
    
    % I assume that all the matrices are defined over the same frequency grid.
    % This function is not immune to errors. Error checks must be
    % implemented by the parent function. The error checks are ommited
    % so that the function can be run on GPU without having to transfer any
    % data to CPU.
    

    % The spectral support (or spatial cutoff) of simulation must be able
    % to accommodate the largest spatial frequency that is captured.
    % Therefore, there is no need to pad either objspecG or Po to find the
    % product of shifted objspecG and Po.

     %Assignment in this manner ensures that new variables are on GPU if Ns
     %is on GPU. Note that Ns is real, whereas objSpectrum and Po are
     %complex.
    img=zeros(size(objSpectrum),'like',Ns); 
    loopvar=ones(1,'like',Ns);
    objspecshifted=zeros(size(objSpectrum),'like',objSpectrum);
    %waith=waitbar(loopvar/Ns,['Computing image intensity by summing over the source:' num2str(100*loopvar/Ns,2) '%']);
    
    
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
    
      % Likewise for nshift
      %--------------------------
      negnshift=single(nshift<0); posnshift=single(nshift>=0);
      nSpectrumStart=negnshift.*(-nshift+1) +...
        posnshift.*1;
      nSpectrumEnd=negnshift.*(L)+ ...
        posnshift.*(L-nshift);
    
      nShiftSpectrumStart=negnshift.*1 +...
        posnshift.*(nshift+1);
      nShiftSpectrumEnd=negnshift.*(L+nshift) +...
        posnshift.*L;
    
    % Need to vectorize the while loop. The idea is:
    % Create a stacked shifted object spectrum with dimensions: (fx,fy,S). 
    % Multiply along (fx,fy) dimensions using Po with bsxfun.
    % Obtain IFFT along (fx,fy) dimensions.
    % Mag. square and sum along S.
    while loopvar<=Ns 
        
        % A parallel computing demo at http://www.mathworks.com/products/demos/parallel-computing/paralleldemo_gpu_mandelbrot/
        % suggests that while loop is more efficient on GPU then for
        % loop.

        %waitbar(loopvar/Ns,waith,['Computing image intensity by summing over the source:' num2str(100*loopvar/Ns,2) '%']);
 
         nidxshifted=colon(nShiftSpectrumStart(loopvar),nShiftSpectrumEnd(loopvar));
         midxshifted=colon(mShiftSpectrumStart(loopvar),mShiftSpectrumEnd(loopvar));
         
         nidxunshifted=colon(nSpectrumStart(loopvar),nSpectrumEnd(loopvar));
         midxunshifted=colon(mSpectrumStart(loopvar),mSpectrumEnd(loopvar));
         
         %
         % objspecshifted=circshift(objSpectrum,[nshift(loopvar) mshift(loopvar)]);
         %
         % circular-shift implemented here or indexed copy implemented in
         % sumOverSource provide the same result and have the same cost.
         % But circshift appears in the profiler when the input to this
         % function are gpuArrays,i.e., it runs on the CPU.
         % Therefore, it is better to use the indexed copy and execute it
         % on GPU.
          objspecshifted(nidxshifted,midxshifted)=objSpectrum(nidxunshifted,midxunshifted);

         imgspectrum=objspecshifted.*Po;
           
           imgpoint=fftshift(ifft2(ifftshift( imgspectrum ))); 
           imgpoint=srcInt(loopvar)*abs(imgpoint).^2;
           img=img+imgpoint;
           
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



% Following is the normalization according to Martin's book. We do not use it, because the
% following leads to divide by zero for dark-field system.
%         normfactor=abs(Po).^2.*abs(IcG); 
%         normfactor=sum(normfactor(:));
%         if(normfactor) %If normafactor == 0, we are imaging with dark-field system.
%             imgG=imgG/sum(normfactor(:));
%         end

end
