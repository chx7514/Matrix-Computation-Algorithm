% xor下的线性方程组求解；lightout游戏的解法
function []=solve_xor(m,n)
    N=m*n;
    T1=diag(ones(m-1,1),1);
    T1=T1+T1'+eye(m);
    T2=diag(ones(n-1,1),1);
    T2=T2+T2'+eye(n);
    A=kron(T2,eye(m))+kron(eye(n),T1)-eye(N);
    b=ones(N,1);
    for i=1:N
        for p=i:N
            if A(p,i)~=0
                break
            end
        end
        A([i p],:)=A([p i],:);
        b([i p])=b([p i]);
        for j=i+1:N
            if A(j,i)
                for k=i:N
                    A(j,k)=xor(A(j,k),A(i,k));
                end
                b(j)=xor(b(j),b(i));
            end
        end
    end
    for i=N:-1:2
        for j=i-1:-1:1
            b(j)=xor((A(j,i)&&b(i)),b(j));
        end
    end
    b=reshape(b,m,n);
    imagesc(b);