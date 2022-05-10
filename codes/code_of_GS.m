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

function [A,R]=CGS2(A)
    [m,n]=size(A);
    R=zeros(n,n);
    R2=zeros(n,n);
    
    R(1,1)=norm(A(1:m,1));
    A(1:n,1)=A(1:n,1)/R(1,1);    
    for i=2:n
        R(1:i-1,i)=A(1:m,1:i-1)'*A(1:m,i);
        A(1:m,i)=A(1:m,i)-A(1:m,1:i-1)*R(1:i-1,i);
        R2(1:i-1,i)=A(1:m,1:i-1)'*A(1:m,i);
        A(1:m,i)=A(1:m,i)-A(1:m,1:i-1)*R2(1:i-1,i);
        R(i,i)=norm(A(1:m,i));
        A(1:m,i)=A(1:m,i)/R(i,i);
        R(1:i-1)=R(1:i-1)+R2(1:i-1);
    end
    
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