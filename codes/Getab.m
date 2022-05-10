function [A,b]=Getab(n)
    m=n-1;
    k=m^2;
    A1=4*eye(m);
    A2=-1*diag(ones(m-1,1),1);
    A1=A1+A2+A2';
    A3=-1*diag(ones(k-m,1),m);
    A=A3+A3';
    for i=0:m-1
        A(m*i+1:m*(i+1),m*i+1:m*(i+1))=A1;
    end
    b=zeros(k,1);
    for i=1:m
        b(i)=sin(i*pi/n);
    end