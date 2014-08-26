R=10;
b=10;
c=1;
pr=.9;
ps=.5;
s=0;
rhoinfo=1.01;
rhoint=1.01;
%%
povals=0:.1:1;
Npo=length(povals);
stabeq=zeros(Npo,3);

for i=1:Npo
    po=povals(i);
    e=equilibria_assortment_interactions(b,c,R,po,rhoinfo,pr,s,ps,rhoint);
    e(2,:)=[];
    e=round(e*10000)/10000;
%     e=[e;1 0 0; 0 1 0; 0 0 1];
    V=stability_coop(b,c,R,po,rhoinfo,pr,s,ps,rhoint);
    V=real(V);
    stable=(sum(V<=0,1)==2);
    if sum(stable)>0
        stabeq(i,:)=e(stable,:);
    end
end
%%
pograd=struct();
pograd.sgn=zeros(Npo,2);
pograd.val=zeros(Npo,2);
for i=1:Npo
    po=povals(i);
    x=stabeq(i,:);
    z=pogradient(x,b,c,R,po,rhoinfo,pr,s,ps,rhoint);
    pograd.val(i,1)=z;
    pograd.sgn(i,1)=sgn(z);
    x=[0 1 0];
    z=pogradient(x,b,c,R,po,rhoinfo,pr,s,ps,rhoint);
    pograd.val(i,2)=z;
    pograd.sgn(i,2)=sgn(z);
end

%%
pip=zeros(Npo,Npo);

for i=1:Npo
    for j=1:Npo
        po=povals(i);
        poprime=povals(j);
        x=stabeq(i,:);
        r = invasionrate_fun(x,b,c,R,po,poprime,rhoinfo,pr,s,ps,rhoint);
        pip(i,j)=r;
    end
end

%%
sign=[-1 1]; 

podelta=.001;
eps=podelta/2;
maxcount=100;
po_init=.9;
postore=zeros(maxcount+1,1);
postore(1)=po_init;
ponow=po_init;
xstore=zeros(maxcount+1,3);
count=0;
while count<=maxcount

    e=equilibria_assortment_interactions(b,c,R,ponow,rhoinfo,pr,s,ps,rhoint);
    x=e(1,:);
    x1=x(1);
    x2=x(2);
    x3=x(3);
    x4=0;
    xstore(count+1,:)=x;
    
    pograd=pogradient(x,b,c,R,ponow,rhoinfo,pr,s,ps,rhoint);
    ponow=ponow+sgn(pograd)*podelta;
    ponow=round(ponow*1000)/1000;
    if ponow>=(0-eps) && ponow<=(1+eps)
        count=count+1;   
        postore(count+1)=ponow;
    else
        postore((count+2):end)=[];
        xstore((count+2):end,:)=[];
        count=maxcount+1;
    end        
end

