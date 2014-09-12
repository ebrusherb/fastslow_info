R=10;
ps=.4;
testi=14;
testj=[1:Npr];
podelta=.00001;
prdelta=.00001;
po=povals(testi);
for j=1:length(testj)
    pr=prvals(testj(j));
    x=stabeq_fun(b,c,R,po,rhoinfo,pr,s,ps,rhoint);
    checkeq(:,j)=col(x);
    checkpayoffs(:,j)=payoffs_assortment(x,b,c,R,po,rhoinfo,pr,s,ps,rhoint);
    
    grad=col(poprgradient(x,b,c,R,po,rhoinfo,pr,s,ps,rhoint));
    checkgrad(:,j)=grad;
    poprime=po;
    prprime=pr+prdelta;
    P = compare_payoffs(x,b,c,R,po,poprime,rhoinfo,pr,prprime,s,ps,rhoint);
    checkdiff(:,j)=P(:,2)-P(:,1);
end

%%
f=find(checkdiff(3,:)>0,1,'first')
for i=1:4
    subplot(2,2,i);
    plot(checkdiff(i,f:end),'-o')
end

%%
ps=.5;
x=[.3 .3 .4];
x1=x(1);
x2=x(2);
x3=x(3);
x4=0;
interaction_matrix_mut([x1 x2 x3 x4],rhoint);
testi=14;
prvals2=0:.01:1;
Npr2=length(prvals2);
testj=[1:(Npr2-2)];
podelta=.1;
prdelta=.01;
po=povals(testi);
checkg=zeros(length(testj),1);
for j=1:length(testj)
    pr=prvals2(testj(j));
    prprime=pr+prdelta;
    g=gbar_assortment_mut(R,po,poprime,rhoinfo,pr,prprime,x1,x3,x4,ps,rhoint);
    checkg(j)=g(1,2)-g(1,1);
end

plot(checkg)