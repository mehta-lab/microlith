function norm=gray2norm(grayscale)
% norm=gray2norm(grayscale) Convert grayscale image to double and normalize
% the maximum value to 1. 
% Useful for comparing contrast in images which have different absolute birghtness.

norm=double(grayscale);
norm=norm./max(abs(norm(:)));
end