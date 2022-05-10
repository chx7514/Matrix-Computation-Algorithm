% Lanczos algorithm for the symmetric eigenvalue problem
function []=Lanczos_for_eigen(A,q,kmax)
    n=size(q,1);
    T=zeros(n);
    Q=zeros(n);
    Q(:,1)=q/norm(q);
    w=A*Q(:,1);
    T(1,1)=Q(:,1)'*w;
    r=w-Q(:,1)*T(1,1);
    T(2,1)=norm(r);
    T(1,2)=T(2,1);
    k=1;
    [~,D]=eig(T(1:k,1:k));
    scatter(k*ones(k,1),diag(D),2,'filled');
    hold on
    while (k<kmax) && T(k+1,k)~=0
        k=k+1;
        Q(:,k)=r/T(k,k-1);
        w=A*Q(:,k);
        T(k,k)=Q(:,k)'*w;
        if k<n
            r=w-Q(:,k)*T(k,k)-Q(:,k-1)*T(k,k-1);
            T(k+1,k)=norm(r);
            T(k,k+1)=T(k+1,k);
        end
        [~,D]=eig(T(1:k,1:k));
        scatter(k*ones(k,1),diag(D),2,'filled');
%         if mod(k,10)==0 %显示正交性
%             subplot(3,4,k/10);
%             imagesc(log10(abs(Q(:,1:k)'*Q(:,1:k)-eye(k))));
%             colorbar;
%             title(['迭代',num2str(k),'次的正交性']);
%         end
    end