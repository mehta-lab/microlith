function etchProfile=mblnnfSiemens(xaxis,yaxis,etchTransition)
% mblnnfSiemens Generate etched profile of MBL/NNF Siemens star.
% All units are in um.
% 
% Author and copyright: 2013. 
%                       Shalin Mehta, www.mshalin.com
%                       Marin Biological Laboratory.
% License: GPL v3.

% Etching extends upto this radius towards the center.
Rin=0.6;
% Outer radius of the etched pattern.
Rout=75/2;

% Number of etched wedges. Equivalently number of wedge pairs.
NWedge=36;

% etchTransition is the distance over which the etched region transitions
% from fully developed to not developed at all. Positive value implies that
% the target was over-exposed and therefore over-etched. Negative value
% implies that the target was under-exposed and therefore under-etched.
% This transition has critical effect on the image contrast measured near
% the center of the pattern, where the transition length can be significant
% fraction or even comparable to the period of the pattern.
% 
% In dark-field towards the center, the wedge where there is no etching
% (traceable from background) is brighter, implying that that part of the
% wedge is sharper than the other. Implying that etched region increases in
% relative width as center is approached, i.e, etching is overefficient.



[xx, yy]=meshgrid(xaxis,yaxis);
[theta,rr]=cart2pol(xx,yy);
dx=xaxis(2)-xaxis(1);
squareProfile=square(NWedge*theta)>0;

etchRadius=round(abs(etchTransition)/dx);

[etchX,etchY]=meshgrid(-etchRadius:etchRadius);
etchR=sqrt(etchX.^2+etchY.^2);
etchNhood=(etchR<etchRadius);
% Flat etch profile: Assume that the pattern can be either etched or not.
etchHeight=double(etchNhood);
% Flat etch profile gives a much better match to experimental dark-field
% image.

% Conical etch profile: Assume that pattern can be partially etched.
% etchHeight=1-(abs(etchR)/etchRadius);
% etchHeight(etchHeight<0)=0;
etchPattern=strel(etchNhood,etchHeight);

switch(sign(etchTransition))
    case -1
        % Assuming that the pattern is under-exposed.
        % The etching is not enough around transition in the mask and
        % therefore unetched part of the pattern is larger.
        squareProfileEroded=imerode(single(squareProfile),etchPattern);
         etchProfile=squareProfileEroded+1; % Correct for offset removed by grayscale erosion.
       % etchProfile=squareProfileEroded;
    case 0
        % ASSUMING THAT THE PATTERN IS EXPOSED JUST RIGHT
        etchProfile=single(squareProfile);

    case 1
        % Assuming that the pattern is over-exposed.
        % The etching is too much around transition in the mask and
        % therefore etched part of pattern is larger.
        squareProfileDilated=imdilate(single(squareProfile),etchPattern);
        etchProfile=squareProfileDilated-1; % Correct for offset added by grayscale dilation.

end

etchProfile=etchProfile.*single(rr<=Rout).*single(rr>=Rin);

end
