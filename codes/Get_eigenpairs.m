%input: a real schur form matrix A
%output: all eigenvalues of A
function [Q,x]=Get_eigenpairs(T,A)
    n=size(T,1);
    x=zeros(n,1);
    i=1;
    Q=zeros(n);
    R=zeros(n,1);
    while 1
        if T(i+1,i)==0
            x(i)=T(i,i);
            Q(:,i)=real_eigen(A,x(i));
            R(i)=norm(A*Q(:,i)-Q(:,i)*x(i));
            i=i+1;
        else 
            b=trace(T(i:i+1,i:i+1));
            c=det(T(i:i+1,i:i+1));
            delta=b*b-4*c;
            x(i)=(b+sqrt(delta))/2;
            Q(:,i)=im_eigen(A,x(i));
            R(i)=norm(A*Q(:,i)-Q(:,i)*x(i));
            R(i+1)=R(i);
            Q(:,i+1)=conj(Q(:,i));
            x(i+1)=(b-sqrt(delta))/2;
            i=i+2;
        end
        if i>=n
            break
        end
    end
    if i==n
        x(n)=T(n,n);
        Q(:,i)=real_eigen(A,x(i));
        R(i)=norm(A*Q(:,i)-Q(:,i)*x(i));
    end
    semilogy(1:n,R);