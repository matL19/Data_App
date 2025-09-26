function [length,pts] = measure_image(image_filepath,image_filename,ref_pts,ref_length,units,varargin)

if numel(varargin) > 0 && numel(varargin) == 2 && varargin{1} == "image scale"
    image_scale = varargin{2};
else
    error("Only valid name/value pair is ""image scale"".")
end

if ~exist('image_scale','var')
    image_scale = 1;
end

cd(image_filepath)
image = imread(image_filename);

figure();
imshow(image*image_scale)
fprintf("You will now select the points to measure. Press enter to continue.\n")
pause
pts = ginput(2);
close
pixels = sqrt(sum((pts(2,:)-pts(1,:)).^2));
[length_per_pixel, fig] = get_distance_per_pixel(image_filepath,image_filename,ref_pts,ref_length,units,'image scale',image_scale);

length = pixels*length_per_pixel;
fprintf("Measured distance: %.10f " + units + "\n",length)
figure(fig.Number)
hold on
plot(pts(:,1),pts(:,2),'co-')


end