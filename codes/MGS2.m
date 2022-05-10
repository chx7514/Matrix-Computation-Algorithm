% 带重正交化的MGS
function [A,R]=MGS2(A)
    [m,n]=size(A);
    R=zeros(n,n);
    R2=zeros(n,n);
    for i=1:n
        for j=1:i-1
            R(j,i)=A(1:m,j)'*A(1:m,i);
            A(1:m,i)=A(1:m,i)-R(j,i)*A(1:m,j);
        end
        for j=1:i-1
            R2(j,i)=A(1:m,j)'*A(1:m,i);
            A(1:m,i)=A(1:m,i)-R2(j,i)*A(1:m,j);
        end
        R(i,i)=norm(A(1:m,i));
        A(1:m,i)=A(1:m,i)/R(i,i);
        R(1:i-1)=R(1:i-1)+R2(1:i-1);
    end