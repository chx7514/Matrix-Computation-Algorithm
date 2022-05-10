% conjugate gradient
function [x]=CG(A,b,x0,kmax)
    x=x0;
    r=b-A*x;
    k=0;
    p=r;
    tol=1e-15;
    norm_r=norm(r);
    res=norm_r;
    while norm_r>norm(b)*tol && k<kmax
        k=k+1;
        t=A*p;
        a=p'*r/(p'*t);
        x=x+p*a;
        l=r'*r;
        r=r-t*a;
        norm_r=norm(r);
        res=[res;norm_r];
        beta=r'*r/l;
        p=r+p*beta;
    end
    figure(1)
    plot([0:k]',res)
    title('convergence history of CG')