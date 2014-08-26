function z = f2(x3) %find equilibrium coop / disc mix

b=evalin('base','b');
c=evalin('base','c');
R=evalin('base','R');
po=evalin('base','po');
rhoinfo=evalin('base','rhoinfo');
pr=evalin('base','pr');
s=evalin('base','s');
ps=evalin('base','ps');
rhoint=evalin('base','rhoint');
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
    