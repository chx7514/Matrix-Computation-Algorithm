% display a point, implemented for the reproduction experiment
function []=displayP(x,y)
    plot(x,y,'black');
    hold on
    scatter(x,y,10,'.r')
    axis([-0.5,0.5,-0.5,0.5]);