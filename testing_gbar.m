x=[.3 .3 .4];
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
X=interaction_matrix_mut([x1 x2 x3 x4],rhoint);
k=zeros(2,4);
G=zeros(2,2);
G2=zeros(2,2);
Gnew=zeros(2,2);
R=4;
if R>1
    for t = 2:R
        k(1,1)=kbar(t-1,pomat(1,1),pr);
        k(1,2)=kbar(t-1,pomat(1,2),pr);
        k(1,3)=kbar(t-1,pomat(1,3),pr);
        k(1,4)=kbar(t-1,pomat(1,4),pr);
        k(2,1)=kbar(t-1,pomat(2,1),prprime);
        k(2,2)=kbar(t-1,pomat(2,2),prprime);
        k(2,3)=kbar(t-1,pomat(2,3),prprime);
        k(2,4)=kbar(t-1,pomat(2,4),prprime);
        Gnew(1,1) = (t-1)*pr*pomat(1,3)*ps + pr*pomat(1,3)*X(3,1)*k(1,1)-pr*pomat(1,3)*ps*(X(3,1)*k(1,1)+X(3,2)*k(1,2)+X(3,3)*k(1,3)+X(3,4)*k(1,4))+pr*pomat(1,3)*X(3,3)*G(1,1)+pr*pomat(1,3)*X(3,4)*G(1,2)+pr*(1-pomat(1,3))*G(1,1);
        Gnew(1,2) = (t-1)*pr*pomat(1,4)*ps + pr*pomat(1,4)*X(4,1)*k(2,1)-pr*pomat(1,4)*ps*(X(4,1)*k(2,1)+X(4,2)*k(2,2)+X(4,3)*k(2,3)+X(4,4)*k(2,4))+pr*pomat(1,4)*X(4,3)*G(2,1)+pr*pomat(1,4)*X(4,4)*G(2,2)+pr*(1-pomat(1,4))*G(1,2);
        Gnew(2,1) = (t-1)*prprime*pomat(2,3)*ps + prprime*pomat(2,3)*X(3,1)*k(1,1)-prprime*pomat(2,3)*ps*(X(3,1)*k(1,1)+X(3,2)*k(1,2)+X(3,3)*k(1,3)+X(3,4)*k(1,4))+prprime*pomat(2,3)*X(3,3)*G(1,1)+prprime*pomat(2,3)*X(3,4)*G(1,2)+prprime*(1-pomat(2,3))*G(2,1);
        Gnew(2,2) = (t-1)*prprime*pomat(2,4)*ps + prprime*pomat(2,4)*X(4,1)*k(2,1)-prprime*pomat(2,4)*ps*(X(4,1)*k(2,1)+X(4,2)*k(2,2)+X(4,3)*k(2,3)+X(4,4)*k(2,4))+prprime*pomat(2,4)*X(4,3)*G(2,1)+prprime*pomat(2,4)*X(4,4)*G(2,2)+prprime*(1-pomat(2,4))*G(2,2);
        G=Gnew;
        
        Gnew2(1,1) = pr*(pomat(1,3)*((t-1)*ps + X(3,1)*k(1,1)-ps*(X(3,1)*k(1,1)+X(3,2)*k(1,2)+X(3,3)*k(1,3)+X(3,4)*k(1,4))+X(3,3)*G2(1,1)+X(3,4)*G2(1,2))+(1-pomat(1,3))*G2(1,1));
        Gnew2(1,2) = pr*(pomat(1,4)*((t-1)*ps + X(4,1)*k(2,1)-ps*(X(4,1)*k(2,1)+X(4,2)*k(2,2)+X(4,3)*k(2,3)+X(4,4)*k(2,4))+X(4,3)*G2(2,1)+X(4,4)*G2(2,2))+(1-pomat(1,4))*G2(1,2));
        Gnew2(2,1) = prprime*(pomat(2,3)*((t-1)*ps + X(3,1)*k(1,1)-ps*(X(3,1)*k(1,1)+X(3,2)*k(1,2)+X(3,3)*k(1,3)+X(3,4)*k(1,4))+X(3,3)*G2(1,1)+X(3,4)*G2(1,2))+(1-pomat(2,3))*G2(2,1));
        Gnew2(2,2) = prprime*(pomat(2,4)*((t-1)*ps + X(4,1)*k(2,1)-ps*(X(4,1)*k(2,1)+X(4,2)*k(2,2)+X(4,3)*k(2,3)+X(4,4)*k(2,4))+X(4,3)*G2(2,1)+X(4,4)*G2(2,2))+(1-pomat(2,4))*G2(2,2));
        G2=Gnew2;
    end
end

t=R+1;
k(1,1)=kbar(t-1,pomat(1,1),pr);
        k(1,2)=kbar(t-1,pomat(1,2),pr);
        k(1,3)=kbar(t-1,pomat(1,3),pr);
        k(1,4)=kbar(t-1,pomat(1,4),pr);
        k(2,1)=kbar(t-1,pomat(2,1),prprime);
        k(2,2)=kbar(t-1,pomat(2,2),prprime);
        k(2,3)=kbar(t-1,pomat(2,3),prprime);
        k(2,4)=kbar(t-1,pomat(2,4),prprime);
        Gfinal(1,1) = (t-1)*pr*pomat(1,3)*ps + pr*pomat(1,3)*X(3,1)*k(1,1)-pr*pomat(1,3)*ps*(X(3,1)*k(1,1)+X(3,2)*k(1,2)+X(3,3)*k(1,3)+X(3,4)*k(1,4))+pr*pomat(1,3)*X(3,3)*G(1,1)+pr*pomat(1,3)*X(3,4)*G(1,2)+pr*(1-pomat(1,3))*G(1,1);
        Gfinal(1,2) = (t-1)*pr*pomat(1,4)*ps + pr*pomat(1,4)*X(4,1)*k(2,1)-pr*pomat(1,4)*ps*(X(4,1)*k(2,1)+X(4,2)*k(2,2)+X(4,3)*k(2,3)+X(4,4)*k(2,4))+pr*pomat(1,4)*X(4,3)*G(2,1)+pr*pomat(1,4)*X(4,4)*G(2,2)+pr*(1-pomat(1,4))*G(1,2);
        Gfinal(2,1) = (t-1)*prprime*pomat(2,3)*ps + prprime*pomat(2,3)*X(3,1)*k(1,1)-prprime*pomat(2,3)*ps*(X(3,1)*k(1,1)+X(3,2)*k(1,2)+X(3,3)*k(1,3)+X(3,4)*k(1,4))+prprime*pomat(2,3)*X(3,3)*G(1,1)+prprime*pomat(2,3)*X(3,4)*G(1,2)+prprime*(1-pomat(2,3))*G(2,1);
        Gfinal(2,2) = (t-1)*prprime*pomat(2,4)*ps + prprime*pomat(2,4)*X(4,1)*k(2,1)-prprime*pomat(2,4)*ps*(X(4,1)*k(2,1)+X(4,2)*k(2,2)+X(4,3)*k(2,3)+X(4,4)*k(2,4))+prprime*pomat(2,4)*X(4,3)*G(2,1)+prprime*pomat(2,4)*X(4,4)*G(2,2)+prprime*(1-pomat(2,4))*G(2,2);
        
        
        Gfinal2(1,1) = pr*(pomat(1,3)*((t-1)*ps + X(3,1)*k(1,1)-ps*(X(3,1)*k(1,1)+X(3,2)*k(1,2)+X(3,3)*k(1,3)+X(3,4)*k(1,4))+X(3,3)*G2(1,1)+X(3,4)*G2(1,2))+(1-pomat(1,3))*G2(1,1));
        Gfinal2(1,2) = pr*(pomat(1,4)*((t-1)*ps + X(4,1)*k(2,1)-ps*(X(4,1)*k(2,1)+X(4,2)*k(2,2)+X(4,3)*k(2,3)+X(4,4)*k(2,4))+X(4,3)*G2(2,1)+X(4,4)*G2(2,2))+(1-pomat(1,4))*G2(1,2));
        Gfinal2(2,1) = prprime*(pomat(2,3)*((t-1)*ps + X(3,1)*k(1,1)-ps*(X(3,1)*k(1,1)+X(3,2)*k(1,2)+X(3,3)*k(1,3)+X(3,4)*k(1,4))+X(3,3)*G2(1,1)+X(3,4)*G2(1,2))+(1-pomat(2,3))*G2(2,1));
        Gfinal2(2,2) = prprime*(pomat(2,4)*((t-1)*ps + X(4,1)*k(2,1)-ps*(X(4,1)*k(2,1)+X(4,2)*k(2,2)+X(4,3)*k(2,3)+X(4,4)*k(2,4))+X(4,3)*G2(2,1)+X(4,4)*G2(2,2)))+prprime*((1-pomat(2,4))*G2(2,2));
        
    