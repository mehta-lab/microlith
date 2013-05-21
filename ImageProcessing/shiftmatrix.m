function out=shiftmatrix(in,xystep,shiftvec,padval)
% OUT = SHIFTMATRIX(IN,XYSTEP,SHIFTVEC,PADVAL) shifts the matrix IN while
% padding the shifted version by PADVAL. XYSTEP specifies the sampling
% interval used to produce the matrix, and should be set to 1 when the
% sampling-interval is not relevant and one wants to shift on pixel-basis.
% SHIFTVEC is a two-element vector ([XSHIFT YSHIFT]) specifying shift along
% the horizontal and vertical, respectively. Positive shift moves the
% matrix to right/down. Ratio of SHIFTVEC to XYSTEP must be integer. If
% PADVAL is set to 'circular', circular shift is employed.  This code is
% roughly 2.5 times faster than using MATLAB's interp2 (with
% nearest-neighbor interpolation) for the task at hand.
%
% Note: The convention of positive shift in Y-direction moving the matrix
% down-wards is OK because, MATLAB's imagesc and imshow plot the matrices
% with  Y-axis running from top to bottom. 
% Author: Shalin Mehta (c) 2010,2011. shalin.mehta@gmail.com


pixshift=shiftvec/xystep; % Amount of shift in pixels.

if( ~all( pixshift-round(pixshift) < 1E-10) ) %Check if all shifts lead to integer pixel-shifts.
    error('The ratio of the amount of shift and the sampling rate should be an integer.');
end

colYshift=round(pixshift(2)); %Due to numerical errors the ratios may not `perfectly' be integer.
rowXshift=round(pixshift(1));

if(strcmp(padval,'circular'))
    out=circshift(in,[colYshift rowXshift]);
else
    out=padval*ones(size(in),class(in)); % Fill the output array with value to be padded.

 % Inspiration for this part comes from snippet on circular shift in: "MATLAB array manipulations tips and
 %tricks" guide by Peter J. Acklam (http://home.online.no/~pjacklam).
 %    inidx = repmat({':'}, ndims(in), 1); % initialize subscripts
 %    outidx= repmat({':'}, ndims(out), 1); 
     coln = size(in, 1); % length along Y-dimension
     rown = size(in, 2); % length along X-dimension.
     
     if(colYshift<0)
         outidx{1}=[1:coln+colYshift];
         inidx{1}=[-colYshift+1:coln];
     else
         outidx{1}=[1+colYshift:coln];
         inidx{1}=[1:coln-colYshift];
     end
    
     if(rowXshift<0)
         outidx{2}=[1:rown+rowXshift];
         inidx{2}=[-rowXshift+1:rown];
     else
         outidx{2}=[1+rowXshift:rown];
         inidx{2}=[1:rown-rowXshift];
     end

     if(size(in,3)>1) %For a stack don't do anything to Z (yet). In future, one may want to shift in Z or in 4th dimension.
         inidx{3}=[1:size(in,3)];
         outidx{3}=inidx{3};
     end
     
     out(outidx{:})=in(inidx{:});
end

end
