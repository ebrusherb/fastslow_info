function z = findx3(x3,b,c,R,po,rhoinfo,pr,s,ps,rhoint)
        x1=1-x3;
        po1=po/(x1+rhoinfo*x3);
        k1=kbar(R,po1,pr);
        X=interaction_matrix([x1 0 x3],rhoint);
        z=b*k1*X(1,3)-c*R;
    end
    