function X = equilibria_assortment_interactions(b,c,R,po,rhoinfo,pr,s,ps,rhoint)
    
    function z = f2(x3) %find equilibrium coop / disc mix
        x1=1-x3;
        po1=po/(x1+rhoinfo*x3);
        po3=po*rhoinfo/(x1+rhoinfo*x3);
        k1=kbar(R,po1,pr);
        k3=kbar(R,po3,pr);
        g=gbar_assortment(R,po,rhoinfo,pr,1-x3,x3,ps,rhoint);
        X=interaction_matrix([x1 0 x3],rhoint);
        p1=b*X(1,1)+b*ps*(R-k1)*X(1,3)+b*k1*X(1,3)-c*R;
        p3=b*X(3,1)+b*ps*(R-k3)*X(3,3)-c*ps*(R-k1*X(3,1)-k3*X(3,3))-c*k1*X(3,1)+(b-c)*g*X(3,3)-s;
        z=p3-p1;
    end
    if sign(f2(0))~=sign(f2(1))
        u=fzero(@f2,[0,1]);
    else u=1;
    end
    X2=[1-u,0,u];
    
    function z = f1(x3) %find equilibrium coop / def mix
        x2=1-x3;
        po2=po/(x2+rhoinfo*x3);
        po3=po*rhoinfo/(x2+rhoinfo*x3);
        k2=kbar(R,po2,pr);
        k3=kbar(R,po3,pr);
        g=gbar_assortment(R,po,rhoinfo,pr,0,x3,ps,rhoint);
        X=interaction_matrix([0 x2 x3],rhoint);
        p2=b*X(2,1)+b*ps*(R-k2)*X(2,3);
        p3=b*X(3,1)+b*ps*(R-k3)*X(3,3)-c*ps*(R-k2*X(3,2)-k3*X(3,3))+(b-c)*g*X(3,3)-s;
        z=p3-p2;
    end
    if sign(f1(0))~=sign(f1(1))
    u=fzero(@f1,[0,1]);
    else u=0;
    end
    X1=[0,1-u,u];
    
    function z = findx3(x3)
        x1=1-x3;
        po1=po/(x1+rhoinfo*x3);
        k1=kbar(R,po1,pr);
        X=interaction_matrix([x1 0 x3],rhoint);
        z=b*k1*X(1,3)-c*R;
    end
    if sign(findx3(0))~=sign(findx3(1))
        x3=fzero(@findx3,[0,1]);
    else x3=1;
    end
    
%     x3=c*R/(b*k1); %find interior equilibrium 
    function z = fi(x1)
        x2=1-x1-x3;
        po1=po/(x1+x2+rhoinfo*x3);
        po2=po/(x1+x2+rhoinfo*x3);
        po3=po*rhoinfo/(x1+x2+rhoinfo*x3);
        k1=kbar(R,po1,pr);
        k2=kbar(R,po2,pr);
        k3=kbar(R,po3,pr);
        g=gbar_assortment(R,po,rhoinfo,pr,x1,x3,ps,rhoint);
        X=interaction_matrix([x1 x2 x3],rhoint);
        p2=b*X(2,1)+b*ps*(R-k2)*X(2,3);
        p3=b*X(3,1)+b*ps*(R-k3)*X(3,3)-c*ps*(R-k1*X(3,1)-k2*X(3,2)-k3*X(3,3))-c*k1*X(3,1)+(b-c)*g*X(3,3)-s;
        z=p3-p2;
    end
    if sign(fi(0))~=sign(fi(1-x3))
    u=fzero(@fi,[0,1-x3]);
    Xi=[u,1-u-x3,x3];
    else Xi=[0 0 1];
    end
    
    X=[X2;X1;Xi]; %return all three equilibria
end
        