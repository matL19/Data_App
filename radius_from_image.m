function [r,dx,other_data] = radius_from_image(filepath,path_to_ref,ref_filename,varargin)
% determines the radius and pinhole displacement from an image
% name-value pairs include:
%   
%   "image scale" - any number
%   
%   "dx definition" - "from center" (default) or "from edge"
%   
%   "radius definition" - "centroid" (default) or "radius of curvature" or
%   "circle fit"
%   
%   "units" - "mm" or "um"


% set default values
scale_factor = 1;
dx_def = "from center";
rad_def = "centroid";
units = "mm";
% read name-value pairs
while numel(varargin) >= 2
    var = varargin{1};
    val = varargin{2};
    switch var
        case "image scale"
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
        otherwise
            error(var + " is an invalid name/value pair keyword.")
    end
    varargin = varargin(3:end);
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
imshow(I)
gel_fig = gcf;
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

% Calculate sample radius

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
    y0 = yC;
    x0 = xC - radius;
elseif rad_def == "circle fit"
    [xC, yC, radius] = circfit(gel_edge(:,1),gel_edge(:,2));
    x0 = xC - radius;
    y0 = yC;
    plot(radius*cos(0:0.01:2*pi)+xC,radius*sin(0:0.01:2*pi)+yC,'Color','red')
end
scatter(xC,yC,'filled','red')
plot([x0 xC],[y0 yC],'green','LineWidth',2)

gcf
gel_center = [xC yC];
radius_pixels = radius;
a = annotation('textbox',[0.1 0.8 0.1 0.1],'String',"r = " + radius_pixels + " pixels");
a.FontSize = 16;
a.Color = 'green';
a.EdgeColor = 'white';

% Calculate pinhole displacement

if dx_def == "from center"
    delta_x = pinhole_center(1) - gel_center(1);
    delta_y = pinhole_center(2) - gel_center(2);
    displacement_pixels = sqrt(delta_x^2 + delta_y^2);
    gcf
    plot([gel_center(1) pinhole_center(1)],[gel_center(2) pinhole_center(2)],'LineWidth',2','Color','cyan')
    b = annotation('textbox',[0.1 0.75 0.1 0.1],'String',"dx = " + displacement_pixels + " pixels");
    b.FontSize = 16;
    b.Color = 'cyan';
    b.EdgeColor = 'white';
elseif dx_def == "from edge"
    distance = zeros(size(gel_edge,1),1);
    for ii = 1:size(gel_edge,1)
        delta_x = pinhole_center(1) - gel_edge(ii,1);
        delta_y = pinhole_center(2) - gel_edge(ii,2);
        distance(ii) = sqrt(delta_x.^2 + delta_y.^2);
    end
    x = [1:numel(distance)]';
    out = struct('x',[],'ydata',[],'yfit',[],'res',[],...
        'fobj',[],'G',[],'O',[]);
   
    tic
    
    %do the fit
    [fobj,G,O] = fit(x,distance,'poly2');
    
    toc
    
    figure(227);clf
    hold on
    plot(x,distance,'bo')
    plot(x,fobj(x),'r-','LineWidth',1.5)
    
    %get results
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
    pinhole_to_edge_fit = out;
    
    [displacement_pixels, min_idx] = min(yfit);
    
    edge_part = gel_edge(min_idx,:);
    
    figure(gel_fig.Number);
    hold on
    plot([edge_part(1) pinhole_center(1)],[edge_part(2) pinhole_center(2)],'LineWidth',2','Color','cyan')
    b = annotation('textbox',[0.1 0.75 0.1 0.1],'String',"dx = " + displacement_pixels + " pixels");
    b.FontSize = 16;
    b.Color = 'cyan';
    b.EdgeColor = 'white';
end

% Convert pixels to length
length_per_pixel = get_distance_per_pixel(path_to_ref,ref_filename,[1476 1939;824 828],2000,"um","image scale",10);
if units == "um"
    dx = displacement_pixels * length_per_pixel;
    r = radius_pixels * length_per_pixel;
elseif units == "mm"
    dx = displacement_pixels * length_per_pixel / 1000;
    r = radius_pixels * length_per_pixel / 1000;
end

a.String = a.String + " " + r + " " + units;
b.String = b.String + " " + dx + " " + units;

% Other parameters output

clear out_struct
out_struct.gel_edge = gel_edge;
out_struct.pinhole_edge = pinhole_edge;
out_struct.gel_center = gel_center;
out_struct.pinhole_center = pinhole_center;
out_struct.dx_definition = dx_def;
out_struct.length_per_pixel = length_per_pixel;
if exist('pinhole_to_edge_fit','var')
    out_struct.dx_fit = pinhole_to_edge_fit;
end
other_data = out_struct;

end