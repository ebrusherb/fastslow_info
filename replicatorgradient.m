function grad = replicatorgradient(x,b,c,R,po,rhoinfo,pr,s,ps,rhoint)
po1=po/(x(1)+x(2)+rhoinfo*x(3));
po2=po/(x(1)+x(2)+rhoinfo*x(3));
po3=po*rhoinfo/(x(1)+x(2)+rhoinfo*x(3));

k1=kbar(R,po1,pr); 
k2=kbar(R,po2,pr);
k3=kbar(R,po3,pr);
g=gbar_assortment(R,po,rhoinfo,pr,x(1),x(3),ps,rhoint);

X=interaction_matrix(x,rhoint);

p=zeros(3,1);
p(1)=b*X(1,1)+b*ps*(R-k1)*X(1,3)+b*k1*X(1,3)-c*R;
p(2)=b*X(2,1)+b*ps*(R-k2)*X(2,3);
p(3)=b*X(3,1)+b*ps*(R-k3)*X(3,3)-c*ps*(R-k1*X(3,1)-k2*X(3,2)-k3*X(3,3))-c*k1*X(3,1)+(b-c)*g*X(3,3)-s;

p=p-p(2);
pbar=sum(p.*x);
advantage=p-pbar;
grad=x.*advantage;
end