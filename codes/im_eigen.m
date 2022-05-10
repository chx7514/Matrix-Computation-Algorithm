function [z]=im_eigen(A,l)
    n=size(A,1);
    p=rand(n,1);
    q=rand(n,1);
    a=real(l);
    b=imag(l);
    T=A-a*eye(n);
    L=T^2+b^2*eye(n);
    k=0;
    while k<2
        r=T*p-b*q;
        x=L\r;
        y=(p-T*x)/b;
        v=complex(x,y);
        z=v/norm(v);
        p=real(z);
        q=imag(z);
        k=k+1;
    end