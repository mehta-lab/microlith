function img=sumOverSource(objSpectrum,Po,Ic)
% SUMOVERSOURCE Simulate partially coherent image by summing over the
% source.
%   IMG=SUMOVERSOURCE(objspec,Po,Ic) computes partially
%   coherent image (IMG) using the spectrum of the specimen (OBJSPECTRUM),
%   pupil of the imaging path (Po), and, intensity distribution of the
%   effective source (Ic). 
%
%   SUMOVERSOURCE uses 'sum over the source' algorithm to compute the
%   partially coherent image, i.e., total image is the sum of the images
%   produced due to individual condenser points. I have placed special
%   emphasis on ensuring accuracy of the computation. The computed images
%   match very well with analytical results (for cases where analyical
%   expressions can be derived) Please refer to two accompanying published
%   test benches : One verifies the accuracy of computation and the other
%   illustrates the usage of the code for computing microscopy images.
%
%   Written by Shalin Mehta, www.mshalin.com 
%   License: GPL v3 or later.


    % I assume that all the matrices are defined over the same frequency grid.
    % This function is not immune to errors. Error checks must be
    % implemented by the parent function. The error checks are ommited
    % so that the function can be run on GPU without having to transfer any
    % data to CPU. 
    objspecshifted=complex(zeros(size(Po),class(Po)));
    img=zeros(size(Po),class(Po));
    
    % The spectral support (or spatial cutoff) of simulation must be able
    % to accommodate the largest spatial frequency that is captured.
    % Therefore, there is no need to pad either objspec or Po to find the
    % product of shifted objspec and Po.
    
    %Find non-zero pixels on condenser and the amount of spectral shift
    %induced by them.
    %------------------------------------------
    [nshift,mshift, SrcInt]=find(Ic);
    

    FreqGridLen=size(Ic,1);
    
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
    
        loopvar=1;
        Ns=numel(nshift);
%        Np=numel(find(Po));
   % msgLen=0;
%        ifftscale=FreqGridLen^2/Np;

    while loopvar<=Ns
        % A parallel computing demo at http://www.mathworks.com/products/demos/parallel-computing/paralleldemo_gpu_mandelbrot/
        % suggests that while loop is more efficient on GPU then for
        % loop.
        
        % Display progress in the console. GUI-based progress bar is
        % avoided so that the function can be compiled as a mex file.
       % progress=unit8(100*loopvar/Ns);
        
        %msg=['Computing image intensity by summing over the source:' int2str(progress) '%'];
        %fprintf(repmat('\b',1,msgLen));
        %fprintf(msg);
        %msgLen=numel(msg);
        
         nidxshifted=colon(nShiftSpectrumStart(loopvar),nShiftSpectrumEnd(loopvar));
         midxshifted=colon(mShiftSpectrumStart(loopvar),mShiftSpectrumEnd(loopvar));
         
         nidxunshifted=colon(nSpectrumStart(loopvar),nSpectrumEnd(loopvar));
         midxunshifted=colon(mSpectrumStart(loopvar),mSpectrumEnd(loopvar));
         
         objspecshifted(nidxshifted,midxshifted)=objSpectrum(nidxunshifted,midxunshifted);
         
         imgspectrum=objspecshifted.*Po;
           
%          imgspectrumsignificant=abs(imgspectrum)>1E-20;
%          ifftscale=FreqGridLen^2/sum(imgspectrumsignificant(:));
          
           %IFFT2 algorithm divides the input by numel(Po). The
           %value at zero index of the output of IFFT2 is equal to the number of nonzero
           %elements,i.e., numel(find(Po)). The above scale compensates
           %for both and ensures that an image of a point produced by a clear
           %circular aprture has a peak value of 1. This normalization allows us to
           %compare images (apple to apple) computed over different grid sizes.
           % Normalizing intensity due to each source point may not be
           % correct. In particular, the term in denominator changes the
           % relative strengths of spatial frequency.
  

%           imgpoint=ifftscale*fftshift(ifft2(ifftshift(imgspectrum ))); 
           imgpoint=fftshift(ifft2(ifftshift(imgspectrum ))); 
           imgpoint=SrcInt(loopvar)*abs(imgpoint).^2;
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
%           
          loopvar=loopvar+1;
        
    end

  %  delete(waith);
% Normalize the image such that the values do not depend on the fineness of
% the source grid.
   
%    img=img/Np;

% Following is the normalization according to Martin's book. It ensures
% that a transparent specimen is imaged with unit intensity.
% normfactor=abs(Po).^2.*abs(Ic); We do not use it, because it leads to
% divide by zero for dark-field system. Instead, through normalizations
% perfomed above, we ensure that image of a point under matched
% illumination is unity. The brightness of all the other configurations is
% relative to this benchmark.
%         


end
