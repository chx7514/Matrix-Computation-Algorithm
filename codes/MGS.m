% the MGS algorimth
function [A,R]=MGS(A)
    [m,n]=size(A);
    R=zeros(n,n);    
    for i=1:n
        R(i,i)=norm(A(1:m,i));
        A(1:m,i)=A(1:m,i)/R(i,i);
        if i<n %更新后面的列
            R(i,i+1:n)=A(1:m,i)'*A(1:m,i+1:n);
            A(1:m,i+1:n)=A(1:m,i+1:n)-A(1:m,i)*R(i,i+1:n);
        end
    end