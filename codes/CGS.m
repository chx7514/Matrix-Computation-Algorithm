% CGS 
function [A,R]=CGS(A)
    [m,n]=size(A);
    R=zeros(n,n);
    R(1,1)=norm(A(1:m,1));
    A(1:n,1)=A(1:n,1)/R(1,1);    
    for i=2:n
        R(1:i-1,i)=A(1:m,1:i-1)'*A(1:m,i);
        A(1:m,i)=A(1:m,i)-A(1:m,1:i-1)*R(1:i-1,i);
        R(i,i)=norm(A(1:m,i));
        A(1:m,i)=A(1:m,i)/R(i,i);
    end