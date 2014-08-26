function grad = pogradient(x,b,c,R,po,rhoinfo,pr,s,ps,rhoint)

sign=[-1 1]; 

podelta=.0001;
tol=1e-10;

x1=x(1);
x2=x(2);
x3=x(3);
x4=0;

pomat=zeros(2,4);
pomat(1,1)=po/(x1+x2+rhoinfo*(x3+x4));
pomat(1,2)=po/(x1+x2+rhoinfo*(x3+x4));
pomat(1,3)=po*rhoinfo/(x1+x2+rhoinfo*(x3+x4));
pomat(1,4)=po*rhoinfo/(x1+x2+rhoinfo*(x3+x4));

invasionrate=zeros(2,1);

for j=1:2
    poprime=po+sign(j)*podelta;

    pomat(2,1)=poprime/(x1+x2+rhoinfo*(x3+x4));
    pomat(2,2)=poprime/(x1+x2+rhoinfo*(x3+x4));
    pomat(2,3)=poprime*rhoinfo/(x1+x2+rhoinfo*(x3+x4));
    pomat(2,4)=poprime*rhoinfo/(x1+x2+rhoinfo*(x3+x4));

    k(1,1)=kbar(R,pomat(1,1),pr);
    k(1,2)=kbar(R,pomat(1,2),pr);
    k(1,3)=kbar(R,pomat(1,3),pr);
    k(1,4)=kbar(R,pomat(1,4),pr);
    k(2,1)=kbar(R,pomat(2,1),pr);
    k(2,2)=kbar(R,pomat(2,2),pr);
    k(2,3)=kbar(R,pomat(2,3),pr);
    k(2,4)=kbar(R,pomat(2,4),pr);

    X=interaction_matrix_mut([x1 x2 x3 x4],rhoint);

    g=gbar_assortment_mut(R,po,poprime,rhoinfo,pr,x1,x3,x4,ps,rhoint);

    p=zeros(4,1);
    p(1)=b*X(1,1)+b*ps*(R-k(1,1))*X(1,3)+b*k(1,1)*X(1,3)+b*ps*(R-k(2,1))*X(1,4)+b*k(2,1)*X(1,4)-c*R;
    p(2)=b*X(2,1)+b*ps*(R-k(1,2))*X(2,3)+b*ps*(R-k(2,2))*X(2,4);
    p(3)=b*X(3,1)+b*ps*(R-k(1,3))*X(3,3)+b*g(1,1)*X(3,3)+b*ps*(R-k(2,3))*X(3,4)+b*g(2,1)*X(3,4)-c*ps*(R-k(1,1)*X(3,1)-k(1,2)*X(3,2)-k(1,3)*X(3,3)-k(1,4)*X(3,4))-c*k(1,1)*X(3,1)-c*g(1,1)*X(3,3)-c*g(1,2)*X(3,4)-s;
    p(4)=b*X(4,1)+b*ps*(R-k(1,4))*X(4,3)+b*g(1,2)*X(4,3)+b*ps*(R-k(2,4))*X(4,4)+b*g(2,2)*X(4,4)-c*ps*(R-k(2,1)*X(4,1)-k(2,2)*X(4,2)-k(2,3)*X(4,3)-k(2,4)*X(4,4))-c*k(2,1)*X(4,1)-c*g(2,1)*X(4,3)-c*g(2,2)*X(4,4)-s;

    invasionrate(j)=p(4)-p(3);

end

if abs(po)<=tol
    grad=invasionrate(2)/podelta;
elseif abs(po-1)<=tol
    grad=-invasionrate(1)/podelta;
else grad=(invasionrate(2)-invasionrate(1))/(2*podelta);
end