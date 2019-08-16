function [BW,maskedImage] = segmentImage(X,segThre,radius)
%Determines number of input arguments, and if <2 (i.e. no radius input),
%sets radius to 2
if (nargin<2) 
    radius=2;
end
% Normalize input data to range in [0,1].
Xmin = min(X(:));
Xmax = max(X(:));
X = (X - Xmin) ./ (Xmax - Xmin);

% Threshold image - manual threshold
if(nargin<2)
segThre=3.568600e-01;
end;
BW = X > segThre;

% Open mask with disk
radius = 2;
decomposition = 0;
se = strel('disk', radius, decomposition);
BW = imopen(BW, se);

% Fill holes
BW = imfill(BW, 'holes');

% Create masked image.
maskedImage = X;
maskedImage(~BW) = 0;