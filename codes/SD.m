% steepest descent
function [x]=SD(A,b,x0)
    r=b-A*x0;
    k=0;
    x=x0;
    tol=1e-10;
    norm_r=norm(r);
    res=norm_r;
    while norm_r>tol 
        k=k+1;
        t=A*r;
        a=r'*r/(r'*t);
        x=x+r*a;
        r=r-t*a;
        norm_r=norm(r);
        res=[res;norm_r];
    end
    figure(1)
    plot([1:k+1]',res)
    title('convergence history of SD')