function s = sgn(x)
divisor=x;
divisor(x==0)=1;
s=abs(x)./divisor;
end