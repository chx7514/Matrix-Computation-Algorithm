% Jacobi diagonalization algorithm for real symmetric matrices
function [A]=Jacobi_Method(A)
    n=size(A,1);
    tol=1e-16;
    d=tol*norm(A,'fro');
    S=A(:);
    k=1;
    cnt=0;
    while off(A)>d
        max_a=0;
        for i=2:n
            for j=1:i-1
                if abs(A(i,j))>max_a
                    p=i;q=j;
                    max_a=abs(A(i,j));
                end
            end
        end %找到最大的A(i,j)
        [c,s]=symSchur2(A,p,q);
        J=Get_J(p,q,c,s,n);
        A=J'*A*J;
        k=k+1;
        S=[S,A(:)];
        cnt=cnt+1;
    end
    figure(1)
    hold on
    x_num=1:k;
    for i=1:n^2
        plot(x_num,S(i,:));
    end
    hold off
    cnt
    
function [sum]=off(A)
    n=size(A,1);
    sum=0;
    for i=2:n
        for j=1:i
            if i~=j
                sum=sum+2*A(i,j)^2;
            end
        end
    end
    sum=sqrt(sum);
    
function [c,s]=symSchur2(A,p,q)
    if A(p,q)~=0
        r=(A(p,p)-A(q,q))/(2*A(p,q));
        if r>=0
            t=1/(r+sqrt(1+r^2));
        else
            t=1/(r-sqrt(1+r^2));
        end
        c=1/sqrt(1+t^2);
        s=c*t;
    end

function [J]=Get_J(p,q,c,s,n)
    J=eye(n);
    J(p,p)=c;
    J(q,q)=c;
    J(p,q)=-s;
    J(q,p)=s;
        