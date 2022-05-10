 % preconditioned CG
function [x]=PCG(A,b,x0,L0,kmax) 
    x=x0;
    k=0;
    n=size(A,1);
    r=b-A*x;
    tol=1e-15;
    norm_r=norm(r);
    res=norm_r;
    while norm_r>norm(b)*tol && k<kmax
        z=r;
        for i=1:n-1
            z(i)=z(i)/L0(i,i);
            z(i+1)=z(i+1)-z(i)*L0(i+1,i);
        end
        z(n)=z(n)/L0(n,n);
        for i=n:-1:2
            z(i)=z(i)/L0(i,i);
            z(i-1)=z(i-1)-z(i)*L0(i,i-1);
        end
        z(1)=z(1)/L0(1,1);
        k=k+1;
        if k==1
            p=z;
            t1=r'*z;
        else
            t2=t1;
            t1=r'*z;
            beta=t1/t2;
            p=z+p*beta;
        end
        w=A*p;
        a=t1/(p'*w);
        x=x+p*a;
        r=r-w*a;
        norm_r=norm(r);
        res=[res;norm_r];
    end
    figure(2)
    plot([0:k]',res)
    title('convergence history of PCG')
end