function out=imtransformAffineMat(in,affinemat,interpolation,varargin)
% out=imtransformAffineMat(in,affinemat,interpolation)
% Wrapper around imtransform that ensures that result of affine transform
% is in the same frame of reference as the input.

% Author and Copyright: Shalin Mehta, HFSP Postdoctoral Fellow
%                           Marine Biological Laboratory, Woods Hole, MA
%                           http://www.mshalin.com
% 

optargs.coordinates='matlab';
optargs.PixSize=1;
if(~isempty(varargin))
    optargs=parsepropval(optargs,varargin{:});
end

switch(optargs.coordinates)
    case 'centered'
        xdata=0.5*[-size(in,2) size(in,2)]*optargs.PixSize;
        ydata=0.5*[-size(in,1) size(in,1)]*optargs.PixSize;
    case 'matlab'
        xdata=[1 size(in,2)]*optargs.PixSize;
        ydata=[1 size(in,1)]*optargs.PixSize;
end

fillval=min(in(:));
out=imtransform(in,maketform('affine',affinemat),interpolation,...
                'XData',xdata,'YData',ydata,'UData',xdata,'VData',ydata,...
                'Size',size(in),'FillValues',fillval);
end
