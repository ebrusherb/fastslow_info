function G=gbar_assortment(R,po,rhoinfo,pr,x1,x3,ps,rhoint) 
    G=0;
    x2=1-x1-x3;
    po1=po/(x1+x2+rhoinfo*x3);
    po2=po/(x1+x2+rhoinfo*x3);
    po3=po*rhoinfo/(x1+x2+rhoinfo*x3);
    X=interaction_matrix([x1 x2 x3],rhoint);

    if R>1
        for t = 2:R
            k1=kbar(t-1,po1,pr);
            k2=kbar(t-1,po2,pr);
            k3=kbar(t-1,po3,pr);
            G = (t-1)*pr*po3*ps + pr*po3*X(3,1)*k1-pr*po3*ps*(X(3,1)*k1+X(3,2)*k2+X(3,3)*k3)+pr*(po3*X(3,3)+1-po3)*G;
        end
    end
    