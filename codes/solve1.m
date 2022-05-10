% 求解 min ||Ax-b||_2^2+\lambda||x||_2^2
function [b]=solve1(A,b,lambda)
    [m,n]=size(A);
    A=[A;sqrt(lambda)*eye(n)];
    b=[b;zeros(n,1)];
    for i=1:n
        for j=1:i            
            [c,s]=Given(A(i,i),A(m+j,i)); %Given变换
            A(i,i)=c*A(i,i)+s*A(m+j,i);
            if i<n
                A(m+j,i+1)=-s*A(i,i+1)+c*A(m+j,i+1); %下一列元素改变
                A(i,i+1)=c*A(i,i+1);                 
            end
            T=[c,s;-s,c]*[b(i);b(m+j)];%处理b
            b(i)=T(1);
            b(m+j)=T(2);
        end
    end
    for i=n:-1:2
        b(i)=b(i)/A(i,i);
        b(i-1)=b(i-1)-A(i-1,i)*b(i);
    end
    b(1)=b(1)/A(1,1);
    b=b(1:n);
    
function [c,s]=Given(a,b)
    if b==0
        c=1;
        s=0;
    else
        if abs(b)>abs(a)
            t=a/b;
            s=1/sqrt(1+t^2);
            c=s*t;
        else
            t=b/a;
            c=1/sqrt(1+t^2);
            s=c*t;
        end
    end
    