function plotGelImageAnalysis(I,r,dx,other_data,units,varargin)

% Regenerates the annotated plot of the gel image using the outputs from
% radius_from_image.m
% 
%     I - the image
% 
%     r - the radius
% 
%     dx - the displacement
% 
%     other_data - the other_data struct output from radius_from_image.m
% 
%     units - the units of r and dx. only accepting mm and um
% 
%     varargin - accepted name/value pairs:
%         
%         "AnnotationPosition" - location of the r = and dx = textbox,
%         normalized to the axes, given as [x y]
% 
%         "FontSize" - font size of the text annotation
%
%         "figure" - figure object
%
%         "axes" - axes object
%

annotation_position = [0 1];
font_size = 16;
gel_fig = figure();
gel_ax = axes(gel_fig);
while numel(varargin) >= 2
    var = varargin{1};
    val = varargin{2};
    switch var
        case "AnnotationPosition"
            annotation_position = val;
        case "FontSize"
            font_size = val;
        case "figure"
            close(gel_fig)
            gel_fig = val;
        case "axes"
            gel_ax = val;
        otherwise
            error("Invalid name/value pair.")
    end
    varargin = varargin(3:end);
end

% unpack the data
radius_pixels = r / other_data.length_per_pixel;
displacement_pixels = dx / other_data.length_per_pixel;
if units == "mm"
    radius_pixels = radius_pixels * 1000;
    displacement_pixels = displacement_pixels * 1000;
end

% display the image
imshow(I)
set(gel_fig,'Position',[813   332   868   615])
hold on

% gel center
scatter(gel_ax,other_data.gel_center(1),other_data.gel_center(2),'filled','red')

% gel edge
scatter(gel_ax,other_data.gel_edge(:,1),other_data.gel_edge(:,2),'filled','green')

% radius line
plot(gel_ax,[other_data.gel_center(1)-radius_pixels other_data.gel_center(1)],[other_data.gel_center(2) other_data.gel_center(2)],'green','LineWidth',2)

% resulting circle
plot(gel_ax,radius_pixels*cos(0:0.01:2*pi)+other_data.gel_center(1),radius_pixels*sin(0:0.01:2*pi)+other_data.gel_center(2),'Color','red')

% pinhole center
scatter(gel_ax,other_data.pinhole_center(1),other_data.pinhole_center(2),'filled','red')

% pinhole edge
scatter(gel_ax,other_data.pinhole_edge(:,1),other_data.pinhole_edge(:,2),'filled','cyan')

% annotate measured radius on plot
a = text(annotation_position(1),annotation_position(2),...
    "r = " + radius_pixels + " pixels",'Units','normalized',...
    'FontSize',font_size,'Color','green','EdgeColor','white',...
    'BackgroundColor','white');

% annotate measured displacement on plot
if other_data.dx_definition == "from center"
    plot(gel_ax,[other_data.gel_center(1) other_data.pinhole_center(1)],[other_data.gel_center(2) other_data.pinhole_center(2)],'LineWidth',2','Color','cyan')
elseif other_data.dx_definition == "from edge"
    
    % find the minimum value of the parabola
    simulated_edge_idxs = 1:0.01:numel(other_data.dx_fit.x);
    displacement_pixels = min(other_data.dx_fit.fobj(simulated_edge_idxs));
    
    % find the point on the edge closest to the pinhole center
    V = [other_data.pinhole_center(1) - other_data.gel_center(1) other_data.pinhole_center(2) - other_data.gel_center(2)];
    unitV = V/norm(V);
    closest_point = other_data.gel_center + unitV * radius_pixels;
    
    % annotate the plot
    plot(gel_ax,[closest_point(1) other_data.pinhole_center(1)],[closest_point(2) other_data.pinhole_center(2)],'LineWidth',2','Color','cyan')
end
b = text(annotation_position(1),annotation_position(2)*0.9,"dx = " + displacement_pixels + " pixels",'Units','normalized',...
    'FontSize',font_size,'Color','cyan','EdgeColor','white',...
    'BackgroundColor','white');

% annotate the plot with lengths
a.String = a.String + " " + r + " " + units;
b.String = b.String + " " + dx + " " + units;


end