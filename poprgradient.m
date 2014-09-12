function grad = poprgradient(x,b,c,R,po,rhoinfo,pr,s,ps,rhoint)

sign=[-1 1]; 

% podelta=.0001;
% prdelta=.0001;
podelta=.01;
prdelta=.01;
tol=1e-10;

invasionrate=zeros(2,2);

for i=1:2
    poprime=po+sign(i)*podelta;
    prprime=pr;
    invasionrate(1,i)=invasionrate_fun(x,b,c,R,po,poprime,rhoinfo,pr,prprime,s,ps,rhoint);
end
for i=1:2
    poprime=po;
    prprime=pr+sign(i)*prdelta;
    invasionrate(2,i)=invasionrate_fun(x,b,c,R,po,poprime,rhoinfo,pr,prprime,s,ps,rhoint);
end

grad=zeros(2,1);
if abs(po)<=tol
    grad(1)=invasionrate(1,2)/podelta;
elseif abs(po-1)<=tol
    grad(1)=-invasionrate(1,1)/podelta;
else grad(1)=(invasionrate(1,2)-invasionrate(1,1))/(2*podelta);
end
if abs(pr)<=tol
    grad(2)=invasionrate(2,2)/prdelta;
elseif abs(pr-1)<=tol
    grad(2)=-invasionrate(2,1)/prdelta;
else grad(2)=(invasionrate(2,2)-invasionrate(2,1))/(2*prdelta);
end