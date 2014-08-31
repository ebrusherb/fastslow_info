function [T,Y] = poprdynamics(tmax,b,c,R,po,rhoinfo,pr,s,ps,rhoint)
info=[po pr];
[T,Y]=ode45(@poprgradient,0:1:tmax,info);

function grad = poprgradient(t,info)
    ponow=info(1);
    prnow=info(2);
    x=stabeq_fun(b,c,R,ponow,rhoinfo,prnow,s,ps,rhoint);
    
    sign=[-1 1]; 

    podelta=.0001;
    prdelta=.0001;
    tol=1e-2;

    invasionrate=zeros(2,2);

    for i=1:2
        poprime=ponow+sign(i)*podelta;
        prprime=prnow;
        invasionrate(1,i)=invasionrate_fun(x,b,c,R,ponow,poprime,rhoinfo,prnow,prprime,s,ps,rhoint);
    end
    for i=1:2
        poprime=ponow;
        prprime=prnow+sign(i)*prdelta;
        invasionrate(2,i)=invasionrate_fun(x,b,c,R,ponow,poprime,rhoinfo,prnow,prprime,s,ps,rhoint);
    end

    grad=zeros(2,1);
    if abs(ponow)<=tol
%         grad(1)=invasionrate(1,2)/podelta;
        grad(1)=max(invasionrate(1,2)/podelta,0);
    elseif abs(ponow-1)<=tol
%         grad(1)=-invasionrate(1,1)/podelta;
        grad(1)=min(-invasionrate(1,1)/podelta,0);
    else grad(1)=(invasionrate(1,2)-invasionrate(1,1))/(2*podelta);
    end
    if abs(prnow)<=tol
%         grad(2)=invasionrate(2,2)/prdelta;
        grad(2)=max(invasionrate(2,2)/prdelta,0);
    elseif abs(prnow-1)<=tol
%         grad(2)=-invasionrate(2,1)/prdelta;
        grad(2)=min(-invasionrate(2,1)/prdelta,0);
    else grad(2)=(invasionrate(2,2)-invasionrate(2,1))/(2*prdelta);
    end

end

end