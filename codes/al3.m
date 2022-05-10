% Reproduce the experiments in “From Random Polygon to Ellipse: An Eigenanalysis” by A. N. Elmachtoub and C. F. Van Loan
function []=al3(n)
t=zeros(n,1);
M=GetMn(n);
for i=1:n-1
    t(i+1)=2*i*pi/n;
end
c=sqrt(2/n)*cos(t);
s=sqrt(2/n)*sin(t);
u=2*pi*rand(1);
v=2*pi*rand(1);
x=cos(u)*c+sin(v)*s;x=x/norm(x,1);
y=cos(v)*c+sin(u)*s;y=y/norm(y,2);
x1=[x;x(1)];
y1=[y;y(1)];
scatter(x1,y1,100,'.r')
title('norm 2 & norm inf');
hold on
for i=1:10
    x=M*x;x=x/norm(x,1);
    y=M*y;y=y/norm(y,2);
    x1=[x;x(1)];
    y1=[y;y(1)];
    if mod(i,2)==0
        scatter(x1,y1,100,'.r')
    else
        scatter(x1,y1,100,'.b')
    end
end
hold off