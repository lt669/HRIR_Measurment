function formatDAFxGraph(surface)

caxis([0 1]);  
view(270,90); 
surface.FaceColor = 'interp'; 
surface.LineStyle = ':';
axis off;
pbaspect([2 1 1]);

end