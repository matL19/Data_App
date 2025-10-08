function [length_per_pixel, reference_fig] = get_distance_per_pixel(image_filepath,image_filename,reference_location,reference_length,units,varargin)
% Finds the length per pixel of an image given a reference object in the
% image.
% image_filepath is the entire path to the image including the filename.
%
% reference_location is a 2x2 matrix where the columns are X,Y coordinates
% of two pixels in the image. The distance between them will be calculated 
% in pixels using the distance formula.
%
% reference_length is the known length between the reference points in the 
% image.
%
% units are the units of your input length. the conversion rate in the
% output will also be in these units.
%
% Variable arguments include:
%   "image scale" - scales up the brightness of the image by multiplying
%   each pixel by the provided value.
%   

while numel(varargin) >= 2
    var = varargin{1};
    val = varargin{2};
    switch var
        case "image scale"
            scale_factor = val;
        otherwise
            error("Invalid name/value pair.")
    end
    varargin = varargin(3:end);     
end
if ~exist('scale_factor','var')
    scale_factor = 1;
end

cd(image_filepath);
I = imread(image_filename);
if numel(size(I)) ~= 2
    I = I(:,:,1:3);
    I = rgb2gray(I);
end
I = I*scale_factor;

if units ~= "um" && units ~= "mm"
    error(units + " is an invalid unit. Currently only ""um"" and ""mm"" are supported.")
end

length_pixels = sqrt(sum((reference_location(:,1) - reference_location(:,2)).^2));

% imshow(I)
% set(gcf,'Position',[1   490   681   457])
% hold on
% plot(reference_location(1,:),reference_location(2,:),'ro-')

length_per_pixel = reference_length/length_pixels;

end