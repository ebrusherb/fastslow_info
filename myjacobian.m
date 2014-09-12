function J=myjacobian(X,b,c,R,po,rhoinfo,pr,s,ps,rhoint)
X=col(X); %NEW
% X=round(X*100000)/100000; %NEW

J=zeros(3,3);
p=payoffs_assortment(X,b,c,R,po,rhoinfo,pr,s,ps,rhoint);
p=p-p(2);
if sum(X==[0;0;1])==3
    J(1,1)=p(1)-p(3);
    J(3,1)=-p(1);
    J(3,3)=-p(3);
end

if sum(X==[1;0;0])==3
    J(1,1)=-p(1);
    J(1,3)=-p(3);
    J(3,3)=p(3)-p(1);
end

if sum(X==[0;1;0])==3
    J(1,1)=p(1);
    J(3,3)=p(3);
end

if sum(X==0)<2
    check=[1;2;3];
    step=.000001;
    for j=[1,3]
       change=(check==j)*step;
       change(2)=-step;
       Xprime=X+change;
       gup=replicatorgradient(Xprime,b,c,R,po,rhoinfo,pr,s,ps,rhoint);
       change=-(check==j)*step;
       change(2)=step;
       Xprime=X+change;
       gdown=replicatorgradient(Xprime,b,c,R,po,rhoinfo,pr,s,ps,rhoint);
       J(:,j)=(gup-gdown)/(2*step);
    end
end
J=J([1,3],[1,3]);
end


        