function [L0,L]=create(n)
    L0=zeros(n);
    L0(1,1)=1;
    for i=2:n
        L0(i,i-1:i)=rand(1,2);
    end
    L=L0;
    for i=3:n
        for j=1:i-2
            if rand(1)<0.1
                L(i,j)=0.2*rand(1)-0.1;
            end
        end
    end
end
        