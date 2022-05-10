function [b]=Getb(x,y,z,n)
    b(1,1)=x;
    for i=2:n-1
        b(i,1)=y;
    end
    b(n)=z;