function [M]=GetMn(n)
    M=zeros(n);
    for i=1:n-1
        M(i,i:i+1)=[1,1];
    end
    M(n,1)=1;
    M(n,n)=1;
    M=0.5*M;