% 高斯消去法（选主元）
function [x]=without(A,b)
    n=size(A,1);
    for i=1:n-1
        A(i+1:n,i)=A(i+1:n,i)/A(i,i);
        A(i+1:n,i+1:n)=A(i+1:n,i+1:n)-A(i+1:n,i)*A(i,i+1:n);
    end
    for j=1:n-1
        b(j+1:n)=b(j+1:n)-b(j)*A(j+1:n,j);
    end
    for j=n:-1:2
        b(j)=b(j)/A(j,j);
        b(1:j-1)=b(1:j-1)-b(j)*A(1:j-1,j);
    end
    b(1)=b(1)/A(1,1);
    x=b;
    