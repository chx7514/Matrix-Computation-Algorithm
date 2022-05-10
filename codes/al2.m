% Reproduce the experiments in “From Random Polygon to Ellipse: An Eigenanalysis” by A. N. Elmachtoub and C. F. Van Loan
function []=al2(n)
    M=GetMn(n);
    x=rand(n,1)-0.5;x=x-sum(x)/n;
    y=rand(n,1)-0.5;y=y-sum(y)/n; 
    x=x/norm(x);
    y=y/norm(y);
    figure(1);
    x1=[x;x(1)];
    y1=[y;y(1)];
    displayP(x1,y1);
    for i=1:200
        x=M*x;
        y=M*y;
        if mod(i,20)==0
            x=x/norm(x);
            y=y/norm(y);
        end
        if i==18
            figure(2);
            x=x/norm(x);
            y=y/norm(y);
            x1=[x;x(1)];
            y1=[y;y(1)];
            displayP(x1,y1);
        end
    end
    figure(3);
    x=x/norm(x);
    y=y/norm(y);
    x1=[x;x(1)];
    y1=[y;y(1)];
    displayP(x1,y1);