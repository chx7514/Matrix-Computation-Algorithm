% Reproduce the experiments in “From Random Polygon to Ellipse: An Eigenanalysis” by A. N. Elmachtoub and C. F. Van Loan
function []=al1(n) 
M=GetMn(n); % generate an n*n matrix
x=rand(n,1)-0.5;
y=rand(n,1)-0.5;
x=x/norm(x);
y=y/norm(y);
figure(1);
x1=[x;x(1)];
y1=[y;y(1)];
displayP(x1,y1);
for i=1:100
    x=M*x;
    y=M*y;
    if i==5||i==20
        if i==5
            figure(2);
        end
        if i==20
            figure(3);
        end
        x1=[x;x(1)];
        y1=[y;y(1)];
        displayP(x1,y1);
    end
end
figure(4);
x1=[x;x(1)];
y1=[y;y(1)];
displayP(x1,y1);