% calculate a Householder transformation
function [v,b]= House (x)
	n=length(x);
    v=zeros(n,1);
    if x==0
        b=0;
        return;
    end
	x=x/norm(x,inf);
	sigma=x(2:n)'*x(2:n);
	v(2:n)=x(2:n);
	if sigma==0
		b=0;
	else
		a=sqrt((x(1))^2+sigma);
		if x(1)<=0
			v(1)=x(1)-a;
		else
			v(1)=-sigma/(x(1)+a);
        end
		b=2*(v(1)^2)/(sigma+(v(1)^2));
		v=v/v(1);        
    end
end