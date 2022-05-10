%  Input: a bidiagonal matrix B, which has a non zero superdiagnoal
%  Output: the result of a Golub-Kahan SVD step on B
function [P,B,Q]=SVD_conver(B)
    n=size(B,1);
    if n==1 % in case B is 1*1
        P=1;
        Q=1;
        return
    end
    for i=1:n-1
        if B(i,i)==0
            % perform a Golub-Kahan SVD step on B(1:i,1:i)
            [P1,B(1:i,1:i),Q1]=SVD_conver(B(1:i,1:i)); 
            % zero the superdiagonal entry in the ith row
            P=eye(n);
            for j=i+1:n
                [c,s]=Givens(B(j,j),B(i,j));
                B(j,j)=c*B(j,j)+s*B(i,j);
                B(i,j)=0;
                Tmp=[P(:,i),P(:,j)]*[c,s;-s,c];
                P(:,i)=Tmp(:,1);
                P(:,j)=Tmp(:,2);
                if j<n
                    Tmp=[c,-s;s,c]*[B(i,j+1);B(j,j+1)];
                    B(i,j+1)=Tmp(1);
                    B(j,j+1)=Tmp(2);
                end
            end
            % perform a Golub-Kahan SVD step on B(i+1:n,i+1:n)
            [P2,B(i+1:n,i+1:n),Q2]=SVD_conver(B(i+1:n,i+1:n));
            P=P*blkdiag(P1,P2);
            Q=blkdiag(Q1,Q2);
            return
        end
    end
    % another case is that B is irreducible
    a=B(n,n)*B(n,n)+B(n-1,n)*B(n-1,n);
    % calculate Wilkinson shift of B'*B
    if n>2
        d=(B(n-1,n-1)*B(n-1,n-1)+B(n-2,n-1)*B(n-2,n-1)-a)/2;
    else
        d=(B(1,1)*B(1,1)-a)/2;
    end
    b=B(n-1,n-1)*B(n-1,n);
    u=a-b*b/(d+sign(d)*sqrt(d*d+b*b));
    y=B(1,1)*B(1,1)-u;
    z=B(1,1)*B(1,2);
    Q=eye(n);
    P=eye(n);
    for k=1:n-1
        [c,s]=Givens(y,z); % Get Q
        if k>1
            B(k-1,k)=c*y+s*z;
        end
        Tmp=B(k:k+1,k:k+1)*[c,-s;s,c];
        B(k:k+1,k+1)=Tmp(:,2);
        y=Tmp(1,1);
        z=Tmp(2,1);
        Q(:,k:k+1)=Q(:,k:k+1)*[c,-s;s,c];
        [c,s]=Givens(y,z); % Get P
        B(k,k)=c*y+s*z;
        if k<n-1
            Tmp=[c,s;-s,c]*B(k:k+1,k+1:k+2);
            B(k+1,k+1:k+2)=Tmp(2,:);
            y=Tmp(1,1);
            z=Tmp(1,2);
        else
            B(k:k+1,k+1)=[c,s;-s,c]*B(k:k+1,k+1);
        end
        P(:,k:k+1)=P(:,k:k+1)*[c,-s;s,c];
    end
end