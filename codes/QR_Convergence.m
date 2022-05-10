%input: an unreduced Hessenberg matrix H
%output: an matrix after one Francis QR step on H, an orthogonal matrix Q that Q'*H*Q=H1
function [H,P]=QR_Convergence(H)
    n=size(H,1);
    if n>2
        m=n-1;
        s=H(m,m)+H(n,n);
        t=H(m,m)*H(n,n)-H(m,n)*H(n,m);
        x=H(1,1)^2+H(1,2)*H(2,1)-s*H(1,1)+t;
        y=H(2,1)*(H(1,1)+H(2,2)-s);
        z=H(2,1)*H(3,2);
        P=eye(n);
        for i=0:n-3
            [v,b]=House([x,y,z]');
            T=b*(v*v');
            P=P*blkdiag(eye(i),eye(3)-T,eye(n-i-3));
            p=max(1,i);
            H(i+1:i+3,p:n)=H(i+1:i+3,p:n)-T*H(i+1:i+3,p:n);
            q=min(n,i+4);
            H(1:q,i+1:i+3)=H(1:q,i+1:i+3)-H(1:q,i+1:i+3)*T;
            x=H(i+2,i+1);
            y=H(i+3,i+1);
            if i<n-3
                z=H(i+4,i+1);
            end
        end
        [v,b]=House([x,y]');
        T=b*(v*v');
        P=P*blkdiag(eye(n-2),eye(2)-T);
        H(n-1:n,n-2:n)=H(n-1:n,n-2:n)-T*H(n-1:n,n-2:n);
        H(1:n,n-1:n)=H(1:n,n-1:n)-H(1:n,n-1:n)*T;
    else
        [v,b]=House(H(:,1));
        T=b*(v*v');
        P=eye(2)-T;
        H=P*H*P;
    end