function [radius_pixels, displacement_pixels, out_struct] = radius_of_curvature_from_image(filepath, varargin)

while numel(varargin) >= 2
    var = varargin{1};
    val = varargin{2};
    switch var
        case "image scale"
            scale_factor = val;
        otherwise
            error(var + " is an invalid name/value pair keyword.")
    end
    varargin = varargin(3:end);
end
if ~exist('scale_factor','var')
    scale_factor = 1;
end

slashes = strfind(filepath,"/");
image_path = convertStringsToChars(filepath);
image_filename = image_path(slashes(end)+1:end);
image_path = image_path(1:slashes(end));
cd(image_path);
I = imread(image_filename);
if numel(size(image)) ~= 2
    I = I(:,:,1:3);
    I = rgb2gray(I);
end
I = I*scale_factor;
figure;
imshow(I)
set(gcf,'Position',[813   332   868   615])
hold on

% ---- When selecting edges, the points you select must be fairly evenly
% spaced. I recommend selecting four points first, one in each cardinal
% direction, then selecting four more points, halfway between each pair of
% cardinal points. These 8 points typically suffice from my testing. ----

% select gel edge
fprintf("Now selecting gel edge...Press enter to continue.\n")
guide = annotation('textbox',[0.2 0 0.1 0.1],'String',"Select gel edge. Press enter to continue.");
guide.FontSize = 18;
guide.Color = 'white';
guide.EdgeColor = 'white';
gel_edge = ginput;
scatter(gel_edge(:,1),gel_edge(:,2),'filled','MarkerFaceColor','green')
% gel_center = [mean(gel_edge(:,1)) mean(gel_edge(:,2))];
% scatter(gel_center(1),gel_center(2),100,'filled','MarkerFaceColor','red')
delete(findall(gcf,'type','annotation'))
% select pinhole edge
fprintf("Now selecting pinhole edge...Press enter to continue.\n")
guide = annotation('textbox',[0.2 0 0.1 0.1],'String',"Select pinhole edge. Press enter to continue.");
guide.FontSize = 18;
guide.Color = 'white';
guide.EdgeColor = 'white';
pinhole_edge = ginput;
scatter(pinhole_edge(:,1),pinhole_edge(:,2),'filled','MarkerFaceColor','cyan')
pinhole_center = [mean(pinhole_edge(:,1)) mean(pinhole_edge(:,2))];
scatter(pinhole_center(1),pinhole_center(2),100,'filled','MarkerFaceColor','red')
delete(findall(gcf,'type','annotation'))
%
% Calculate values
delete(findall(gcf,'type','annotation'))
% calculate the radius of curvature
x = gel_edge(:,1);
y = gel_edge(:,2);
x = x';
y = y';
x_dot = gradient(x);
y_dot = gradient(y);
x_double_dot = gradient(x_dot);
y_double_dot = gradient(y_dot);

% formula from wikipedia
radius = abs((x_dot.^2+y_dot.^2).^(3/2)/(x_dot.*y_double_dot-y_dot.*x_double_dot));

if mod(numel(x),2) == 0
    x0 = (x(numel(x)/2)+x(numel(x)/2+1))/2;
    y0 = (y(numel(y)/2)+y(numel(y)/2+1))/2;
else
    x0 = x(numel(x)/2+0.5);
    y0 = y(numel(y)/2+0.5);
end
m = (y(end)-y(1))/(x(end)-x(1));
delta_x = abs(radius/sqrt(1+m^(-2)));
delta_y = abs(-delta_x/m);
chord = [x(end) x(1); y(end) y(1)];
chord_midpt = mean(chord,2);
if x0 > chord_midpt(1)
    xC = x0-delta_x;
else
    xC = x0+delta_x;
end
if y0 > chord_midpt(2)
    yC = y0-delta_y;
else
    yC = y0+delta_y;
end
plot(chord(1,:),chord(2,:),'red')
scatter(x0,y0,'filled','red')
scatter(xC,yC,'filled','red')
plot([x0 xC],[y0 yC],'green','LineWidth',2)

gcf
gel_center = [xC yC];
radius_pixels = radius;
a = annotation('textbox',[0.1 0.8 0.1 0.1],'String',"r = " + radius_pixels + " pixels");
a.FontSize = 16;
a.Color = 'green';
a.EdgeColor = 'white';

% find pinhole displcement
delta_x = pinhole_center(1) - gel_center(1);
delta_y = pinhole_center(2) - gel_center(2);
displacement_pixels = sqrt(delta_x^2 + delta_y^2);
gcf
plot([gel_center(1) pinhole_center(1)],[gel_center(2) pinhole_center(2)],'LineWidth',2','Color','cyan')
b = annotation('textbox',[0.1 0.75 0.1 0.1],'String',"dx = " + displacement_pixels + " pixels");
b.FontSize = 16;
b.Color = 'cyan';
b.EdgeColor = 'white';

% if units == "mm"
%     sub_image = I(scale_bar_region(1,1):scale_bar_region(1,2),scale_bar_region(2,1):scale_bar_region(2,2));
%     threshold = 180;
%     scale_bar = sub_image < threshold;
%     scale_bar = mean(scale_bar,1);
%     scale_bar_pixels = numel(scale_bar(scale_bar > 0));
%     
%     mm_per_pixel = scale_bar_mm/scale_bar_pixels;
%     gcf
%     radius_mm = radius_pixels * mm_per_pixel;
%     a.String = a.String + ", " + radius_mm + " mm";
%     displacement_mm = displacement_pixels * mm_per_pixel;
%     b.String = b.String + ", " + displacement_mm + " mm";
% else
%     radius_mm = radius_pixels;
%     displacement_mm = displacement_pixels;
% end

% other data output structure
clear out_struct
out_struct.gel_edge = gel_edge;
out_struct.pinhole_edge = pinhole_edge;
out_struct.gel_center = gel_center;
out_struct.pinhole_center = pinhole_center;

end