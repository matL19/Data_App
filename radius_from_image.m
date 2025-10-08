function [r,dx,other_data] = radius_from_image(filepath,path_to_ref,ref_filename,varargin)


% [r,dx,other_data_struct] = radius_from_image(path,ref_path,ref_file,varargin)
% 
% determines the radius and pinhole displacement from an image
% name-value pairs include:
% 
%   "image scale" - any integer or float > 0
% 
%   "dx definition" - "from center" (default) or "from edge"
% 
%   "radius definition" - "centroid" (default) or "radius of curvature" or
%   "circle fit"
% 
%   "units" - "mm" or "um"
% 
%   "flag_plot" - true or false



% set default values
scale_factor = 1;
dx_def = "from center";
rad_def = "centroid";
units = "mm";
flag_plot = true;

% read name-value pairs
while numel(varargin) >= 2
    var = varargin{1};
    val = varargin{2};
    switch var
        case "image scale"
            if val <= 0
                error("Value for image scale must be a positive number.")
            end
            scale_factor = val;
        case "dx definition"
            if val == "from center" || val == "from edge"
                dx_def = val;
            else
                error("Only value arguments ""from center"" and ""from edge"" can pair with ""dx definition""");
            end
        case "radius definition"
            if val == "centroid" || val == "radius of curvature" || val == "circle fit"
                rad_def = val;
            else
                error("Only value arguments ""centroid"" and ""radius of curvature"" can pair with ""radius definition""");
            end
        case "units"
            if val == "um" || val == "mm"
                units = val;
            else
                erro("Only um and mm units supported.")
            end
        case "flag_plot"
            flag_plot = val;
        otherwise
            error(var + " is an invalid name/value pair keyword.")
    end
    varargin = varargin(3:end);
end

% get the parts of the filepath
[image_path, image_filename, ext] = fileparts(filepath);

% read the image
cd(image_path);
I = imread(image_filename + ext);

% greyscale the image if in color
if numel(size(I)) ~= 2
    I = I(:,:,1:3);
    I = rgb2gray(I);
end

% display the image for user input
I = I*scale_factor;
gel_fig = figure();
imshow(I)
set(gel_fig,'Position',[813   332   868   615])
hold on

% ---- When selecting edges, the points you select must be fairly evenly
% spaced. I recommend selecting four points first, one in each cardinal
% direction, then selecting four more points, halfway between each pair of
% cardinal points. These 8 points typically suffice from my testing. ----

% select gel edge
fprintf("Now selecting gel edge...Press enter to continue.\n")
guide = annotation('textbox',[0.2 0 0.1 0.1],'String',...
    sprintf("Select gel edge. Current method: %s. Press enter to continue.",rad_def),...
    'BackgroundColor','white','Color','black','FontSize',18,'EdgeColor','none');
gel_edge = ginput;
% delete(findall(gcf,'type','annotation'))

% select pinhole edge
fprintf("Now selecting pinhole edge...Press enter to continue.\n")
guide = annotation('textbox',[0.2 0 0.1 0.1],'String',...
    "Select pinhole edge. Current method: centroid. Press enter to continue.",...
    'BackgroundColor','white','Color','black','FontSize',18,'EdgeColor','none');
pinhole_edge = ginput;
% delete(findall(gcf,'type','annotation'))

close(gel_fig.Number)

% Calculate sample radius
[radius_pixels, gel_center] = calculateRadius(gel_edge, rad_def);

% Calculate pinhole displacement
% produces dx value "displacement_pixels"
[displacement_pixels, pinhole_center, pinhole_to_edge_fit] = calculateDisplacement(pinhole_edge, gel_edge, gel_center, dx_def);

% Convert pixels to length
microns_per_pixel = get_distance_per_pixel(path_to_ref,ref_filename,[1476 1939;824 828],2000,"um","image scale",10);

dx_um = displacement_pixels * microns_per_pixel;
r_um = radius_pixels * microns_per_pixel;
if units == "mm"
    dx = dx_um / 1000;
    r = r_um / 1000;
else
    dx = dx_um;
    r = r_um;
end

% ALL other parameters output
clear other_data
other_data.gel_edge = gel_edge;
other_data.pinhole_edge = pinhole_edge;
other_data.gel_center = gel_center;
other_data.pinhole_center = pinhole_center;
other_data.dx_definition = dx_def;
other_data.radius_definition = rad_def;
other_data.length_per_pixel = microns_per_pixel;
other_data.units = units;
other_data.dx_fit = pinhole_to_edge_fit;

if flag_plot
    
    if other_data.units == "mm"
        plotGelImageAnalysis(I,r,dx,other_data,units)
    elseif other_data.units == "um"
        plotGelImageAnalysis(I,r_um,dx_um,other_data,units)
    end
    
end

end

function [radius_pixels, gel_center] = calculateRadius(gel_edge,rad_def)

if rad_def == "radius of curvature"
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
elseif rad_def == "centroid"
    xC = mean(gel_edge(:,1));
    yC = mean(gel_edge(:,2));
    distances = zeros(1,size(gel_edge,1));
    for ii = 1:size(gel_edge,1)
        delta_x = xC - gel_edge(ii,1);
        delta_y = yC - gel_edge(ii,2);
        distances(ii) = sqrt(delta_x.^2 + delta_y.^2);
    end
    radius = mean(distances);
elseif rad_def == "circle fit"
    [xC, yC, radius] = circfit(gel_edge(:,1),gel_edge(:,2));
end

radius_pixels = radius;
gel_center = [xC yC];

end

function [displacement_pixels, pinhole_center, pinhole_to_edge_fit] = calculateDisplacement(pinhole_edge, gel_edge, gel_center, dx_def)

% calculate the pinhole center using centroid
pinhole_center = [mean(pinhole_edge(:,1)) mean(pinhole_edge(:,2))];

if dx_def == "from center"
    
    delta_x = pinhole_center(1) - gel_center(1);
    delta_y = pinhole_center(2) - gel_center(2);
    
    displacement_pixels = sqrt(delta_x^2 + delta_y^2);

    pinhole_to_edge_fit = struct();
    
elseif dx_def == "from edge"
    
    % calculate the distance from each gel edge point to pinhole center
    distance = zeros(size(gel_edge,1),1);
    for ii = 1:size(gel_edge,1)
        delta_x = pinhole_center(1) - gel_edge(ii,1);
        delta_y = pinhole_center(2) - gel_edge(ii,2);
        distance(ii) = sqrt(delta_x.^2 + delta_y.^2);
    end
    
    % fit this to a parabola
    x = [1:numel(distance)]';
    out = struct('x',[],'ydata',[],'yfit',[],'res',[],...
        'fobj',[],'G',[],'O',[]);
    
    [fobj,G,O] = fit(x,distance,'poly2');
    
    % get fit results
    yfit = fobj(x);
    out.x = 1:numel(distance);
    out.ydata = distance;
    out.yfit = yfit;
    out.res = distance - yfit;
    out.fobj = fobj;
    out.G = G;
    out.O = O;
    
    if out.O.exitflag < 1
        warning('Curve fit did not converge!!! Results might not be trustworthy.');
    end
    
    displacement_pixels = min(out.fobj(1:0.01:numel(out.x)));
    
    pinhole_to_edge_fit = out;
    
  
end

end