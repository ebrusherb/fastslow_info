function r = invasionrate_fun(x,b,c,R,po,poprime,rhoinfo,pr,prprime,s,ps,rhoint)
x1=x(1);
x2=x(2);
x3=x(3);
x4=0;

pomat=zeros(2,4);
pomat(1,1)=po/(x1+x2+rhoinfo*(x3+x4));
pomat(1,2)=po/(x1+x2+rhoinfo*(x3+x4));
pomat(1,3)=po*rhoinfo/(x1+x2+rhoinfo*(x3+x4));
pomat(1,4)=po*rhoinfo/(x1+x2+rhoinfo*(x3+x4));

pomat(2,1)=poprime/(x1+x2+rhoinfo*(x3+x4));
pomat(2,2)=poprime/(x1+x2+rhoinfo*(x3+x4));
pomat(2,3)=poprime*rhoinfo/(x1+x2+rhoinfo*(x3+x4));
pomat(2,4)=poprime*rhoinfo/(x1+x2+rhoinfo*(x3+x4));

k(1,1)=kbar(R,pomat(1,1),pr);
k(1,2)=kbar(R,pomat(1,2),pr);
k(1,3)=kbar(R,pomat(1,3),pr);
k(1,4)=kbar(R,pomat(1,4),pr);
k(2,1)=kbar(R,pomat(2,1),prprime);
k(2,2)=kbar(R,pomat(2,2),prprime);
k(2,3)=kbar(R,pomat(2,3),prprime);
k(2,4)=kbar(R,pomat(2,4),prprime);

% cost=s;
% costprime=s;
totalk=sum(k.*[x1 x2 x3 x4;x1 x2 x3 x4],2);
cost=s*totalk(1);
costprime=s*totalk(2);

X=interaction_matrix_mut([x1 x2 x3 x4],rhoint);

g=gbar_assortment_mut(R,po,poprime,rhoinfo,pr,prprime,x1,x3,x4,ps,rhoint);

p=zeros(4,1);
p(1)=b*X(1,1)+b*ps*(R-k(1,1))*X(1,3)+b*k(1,1)*X(1,3)+b*ps*(R-k(2,1))*X(1,4)+b*k(2,1)*X(1,4)-c*R;
p(2)=b*X(2,1)+b*ps*(R-k(1,2))*X(2,3)+b*ps*(R-k(2,2))*X(2,4);
p(3)=b*X(3,1)+b*ps*(R-k(1,3))*X(3,3)+b*g(1,1)*X(3,3)+b*ps*(R-k(2,3))*X(3,4)+b*g(2,1)*X(3,4)-c*ps*(R-k(1,1)*X(3,1)-k(1,2)*X(3,2)-k(1,3)*X(3,3)-k(1,4)*X(3,4))-c*k(1,1)*X(3,1)-c*g(1,1)*X(3,3)-c*g(1,2)*X(3,4)-cost;
p(4)=b*X(4,1)+b*ps*(R-k(1,4))*X(4,3)+b*g(1,2)*X(4,3)+b*ps*(R-k(2,4))*X(4,4)+b*g(2,2)*X(4,4)-c*ps*(R-k(2,1)*X(4,1)-k(2,2)*X(4,2)-k(2,3)*X(4,3)-k(2,4)*X(4,4))-c*k(2,1)*X(4,1)-c*g(2,1)*X(4,3)-c*g(2,2)*X(4,4)-costprime;

r=p(4)-p(3);
