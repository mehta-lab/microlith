function h=imagecat(varargin)
% ImageCat displays a 'catalog' of images on individual axes and optionally
% links the axes.
%
% IMAGECAT(I1,...) generates a simultaneous display of multiple images on
% individual axes.
% 
% IMAGECAT(x,y,I1,...) displays all the images over axes defined by vectors
% x and y.
% 
% IMAGECAT(...,'link') links all image axes (see linkaxes).
% IMAGECAT(...,'xy') displays all axes with xy orientation (see axis xy)
% IMAGECAT(...,'equal') sets the axes to be proportional (see axis equal)
%
% First two arguments (x,y) and last three arguments ('link','xy','equal')
% are optional. Last three arguments can be issued in any order after the
% image arguments.
% 
% Usage example:
% 
% See also:
% montage, axis, linkaxes.
%
% Author: Shalin Mehta (www.mshalin.com) 
% License: BSD 
% Version history: March 2011, initial implementation.
%                  June 2012, optional linking, axis xy, axis equal.
%                  Feb 2013, improved error checking and included usage example.

clf;

if(nargin < 2)
    error('Atleast two images are required.');
end

% The first two arguments may be axes of the images.
% imageidx lists the indices of images in varargin.
if(isvector(varargin{1}) && isvector(varargin{2}))
    x=varargin{1};
    y=varargin{2};
    imageidx=3:length(varargin); 
else
    x=1:size(varargin{1},2);
    y=1:size(varargin{1},1);
    imageidx=1:length(varargin);
end


% The last four arguments may be 'link', 'xy', 'equal','colormap'.
postargsStart=nargin-4;
if(postargsStart<3) % First two arguments are images or vectors.
    postargsStart=3;
end
    
if(any(cellfun(@(x) strcmp(x,'link'), varargin(postargsStart:end))))
    linkplots=true;
    imageidx=imageidx(1:end-1);
else 
    linkplots=false;
end

if(any(cellfun(@(x) strcmp(x,'xy'), varargin(postargsStart:end))))
    plotxy=true;
    imageidx=imageidx(1:end-1);
else
    plotxy=false;
end

if(any(cellfun(@(x) strcmp(x,'equal'), varargin(postargsStart:end))))
    axiseq=true;
    imageidx=imageidx(1:end-1);
else
    axiseq=false;
end

if(any(cellfun(@(x) strcmp(x,'colorbar'), varargin(postargsStart+1:end))))
    showcolormap=true;
    imageidx=imageidx(1:end-1);
else
    showcolormap=false;
end

% Make sure that all other arguments are images.
ismatrixvector=cellfun(@(x) ndims(x)<=3 & ( size(x,3) == 1 | size(x,3) == 3),...
    varargin(imageidx));
if(~all(ismatrixvector))
    error('One or more of the inputs that should be image are not images.');
end

h=zeros(size(imageidx));
numcols=ceil(length(imageidx)/2);

for n=1:length(imageidx)
    h(n)=subplot(2,numcols,n); 
    
    imagesc(x,y,varargin{imageidx(n)}); 
    title(inputname(imageidx(n)));
    
    if(showcolormap)
        colorbar;
    end
    
    if(plotxy)
        axis xy;
    end
    
    if(axiseq);
        axis equal;
    end
    
    axis tight;
end

if(linkplots)
    linkaxes(h);
end

end