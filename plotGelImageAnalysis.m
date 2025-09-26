function plotGelImageAnalysis(I,r,dx,other_data,units)

% must give r and dx in the same units that are indicated in the struct

% unpack the data
radius_pixels = r / other_data.length_per_pixel;
displacement_pixels = dx / other_data.length_per_pixel;
if units == "mm"
    radius_pixels = radius_pixels * 1000;
    displacement_pixels = displacement_pixels * 1000;
end

% display the image
imshow(I)
gel_fig = gcf;
set(gcf,'Position',[813   332   868   615])
hold on

% gel center
scatter(other_data.gel_center(1),other_data.gel_center(2),'filled','red')

% gel edge
scatter(other_data.gel_edge(:,1),other_data.gel_edge(:,2),'filled','green')

% radius line
plot([other_data.gel_center(1)-radius_pixels other_data.gel_center(1)],[other_data.gel_center(2) other_data.gel_center(2)],'green','LineWidth',2)

% resulting circle
plot(radius_pixels*cos(0:0.01:2*pi)+other_data.gel_center(1),radius_pixels*sin(0:0.01:2*pi)+other_data.gel_center(2),'Color','red')

% pinhole center
scatter(other_data.pinhole_center(1),other_data.pinhole_center(2),'filled','red')

% pinhole edge
scatter(other_data.pinhole_edge(:,1),other_data.pinhole_edge(:,2),'filled','cyan')

% annotate measured radius on plot
a = annotation('textbox',[0.1 0.8 0.1 0.1],'String',"r = " + radius_pixels + " pixels");
a.FontSize = 16;
a.Color = 'green';
a.EdgeColor = 'white';

% annotate measured displacement on plot
if other_data.dx_definition == "from center"
    plot([other_data.gel_center(1) other_data.pinhole_center(1)],[other_data.gel_center(2) other_data.pinhole_center(2)],'LineWidth',2','Color','cyan')
    b = annotation('textbox',[0.1 0.75 0.1 0.1],'String',"dx = " + displacement_pixels + " pixels");
    b.FontSize = 16;
    b.Color = 'cyan';
    b.EdgeColor = 'white';
elseif other_data.dx_definition == "from edge"
    
    % find the minimum value of the parabola
    simulated_edge_idxs = 1:0.01:numel(other_data.dx_fit.x);
    displacement_pixels = min(other_data.dx_fit.fobj(simulated_edge_idxs));
    
    % find the point on the edge closest to the pinhole center
    V = [other_data.pinhole_center(1) - other_data.gel_center(1) other_data.pinhole_center(2) - other_data.gel_center(2)];
    unitV = V/norm(V);
    closest_point = other_data.gel_center + unitV * radius_pixels;
    
    % annotate the plot
    figure(gel_fig.Number);
    hold on
    plot([closest_point(1) other_data.pinhole_center(1)],[closest_point(2) other_data.pinhole_center(2)],'LineWidth',2','Color','cyan')
    b = annotation('textbox',[0.1 0.75 0.1 0.1],'String',"dx = " + displacement_pixels + " pixels");
    b.FontSize = 16;
    b.Color = 'cyan';
    b.EdgeColor = 'white';
end

% annotate the plot with lengths
a.String = a.String + " " + r + " " + units;
b.String = b.String + " " + dx + " " + units;


end