% 高斯消去法（选主元）
function [x]=with(A,b)
    n=size(A,1);
    u=zeros(n-1);
    for i=1:n-1
        max=abs(A(i,i));
        p=i;
        for j=i+1:n
            if abs(A(j,i))>max
                max=abs(A(j,i));
                p=j;
            end
        end
        u(i)=p;
        A([i p],:)=A([p i],:);
        b([i p])=b([p i]);
        if A(i,i)~=0
            A(i+1:n,i)=A(i+1:n,i)/A(i,i);
            A(i+1:n,i+1:n)=A(i+1:n,i+1:n)-A(i+1:n,i)*A(i,i+1:n);
        else
            stop
        end
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
    for i=1:n-1
        x([i u(i)])=x([u(i) i]);
    end