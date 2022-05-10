function [A,d]=HouseQR(A)
    [m,n]=size(A);
    d=zeros(min(m-1,n),1);   
    for j=1:n
        if j<m
            [v,b]=House(A(j:m,j));
            A(j:m,j:n)=A(j:m,j:n)-b*(v*v')*A(j:m,j:n);
            d(j)=b;
            A(j+1:m,j)=v(2:m-j+1);
        end
    end