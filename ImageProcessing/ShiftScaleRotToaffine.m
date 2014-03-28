function [affinemat] = ShiftScaleRotToaffine( xshift,yshift, xscale,yscale,rotation )
% It is assumed that co-efficients of the affine matrix represent the
% transformation in the order: Rotation, Scaling, Translation.
% That is the affine matrix is assumed to be R*S*T
% R.S.T=
% ( sxCost	-sySint	0
%   sxSint	syCost	0
%   xt      yt      1)

affinemat=[xscale*cos(rotation) -yscale*sin(rotation)   0;...
           xscale*sin(rotation)  yscale*cos(rotation)   0;...
           xshift                yshift                 1];

% Little slower version of above is.
% T=@(xt,yt) [1 0 0;0 1 0; xt yt 1];
% R=@(theta) [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
% S=@(sx,sy) [sx 0 0; 0 sy 0; 0 0 1];
% affinemat=R(rotation)*S(xscale,yscale)*T(xshift,yshift);
end

