R=10;
b=10;
c=1;
ps=.5;
s=0.00;
rhoinfo=1.01;
rhoint=1.01;
%with rhoint=0 why is stability not continuous wrt po?!?
%%
povals=0:.1:1;
Npo=length(povals);
prvals=0:.1:1;
Npr=length(prvals);
stabeq=zeros(Npo,Npr,3);

for i=1:Npo
    for j=1:Npr
    po=povals(i);
    pr=prvals(j);
    e=equilibria_assortment_interactions(b,c,R,po,rhoinfo,pr,s,ps,rhoint);
    e(2,:)=[];
    e=round(e*10000)/10000;
%     e=[e;1 0 0; 0 1 0; 0 0 1];
    V=stability_coop(b,c,R,po,rhoinfo,pr,s,ps,rhoint);
    V=real(V);
    stable=(sum(V<=0,1)==2);
    if sum(stable)>0
        stabeq(i,j,:)=e(find(stable>0,1,'first'),:);
    end
    end
end

%%
grad2d=zeros(Npo,Npr,2);

for i=1:Npo
    for j=1:Npr
        po=povals(i);
        pr=prvals(j);
        x=stabeq(i,j,:);
        grad = poprgradient(x,b,c,R,po,rhoinfo,pr,s,ps,rhoint);
        grad2d(j,i,:)=grad;
    end
end
[xcoord,ycoord]=meshgrid(povals,prvals);
xcoord=col(xcoord);
ycoord=col(ycoord);
grad2d_nanrm=reshape(grad2d,[],2);
grad2d_nanrm(isnan(col(grad2d(:,:,1))),:)=[];
xcoord(isnan(col(grad2d(:,:,1))))=[];
ycoord(isnan(col(grad2d(:,:,1))))=[];
%%
figure
scale=1;
imagesc(povals,prvals,transpose(stabeq(:,:,1)))
set(gca,'YDir','normal')
set(gca,'xlim',[0 1],'ylim',[0 1])
hold on
quiver(xcoord,ycoord,scale*sgn(grad2d_nanrm(:,1)),scale*sgn(grad2d_nanrm(:,2)),'Color','white')
hold off

%% mutual invasibility ?!? looks like not.
mutinvasion=[];
t=col(1:(Npo*Npr)^2);
t=t(find(col(tril(ones(Npo*Npr,Npo*Npr),-1)==1)));
Nt=length(t);

for c=1:Nt
    [ind1, ind2]=ind2sub([Npo*Npr, Npo*Npr],t(c));
    [i,j]=ind2sub([Npo,Npr],ind1);
    [k,l]=ind2sub([Npo,Npr],ind2);
    po=povals(i);
    pr=prvals(j);
    x=stabeq(i,j,:);
    poprime=povals(k);
    prprime=prvals(l);
    r1=invasionrate_fun(x,b,c,R,po,poprime,rhoinfo,pr,prprime,s,ps,rhoint);
    x=stabeq(k,l,:);
    r2=invasionrate_fun(x,b,c,R,poprime,po,rhoinfo,prprime,pr,s,ps,rhoint);
    if r1>0 && r2>0
        mutinvasion=[mutinvasion;[po pr poprime prprime]];
    end
end

                
%%
pograd=struct();
pograd.sgn=zeros(Npo,2);
pograd.val=zeros(Npo,2);
for i=1:Npo
    po=povals(i);
    pr=prvals(end-1);
    x=stabeq(i,end-1,:);
    z=pogradient(x,b,c,R,po,rhoinfo,pr,s,ps,rhoint);
    pograd.val(i,1)=z;
    pograd.sgn(i,1)=sgn(z);
    x=[0 1 0];
    z=pogradient(x,b,c,R,po,rhoinfo,pr,s,ps,rhoint);
    pograd.val(i,2)=z;
    pograd.sgn(i,2)=sgn(z);
end

% %%
% pip=zeros(Npo,Npo);
% 
% for i=1:Npo
%     for j=1:Npo
%         po=povals(i);
%         poprime=povals(j);
%         x=stabeq(i,:);
%         r = invasionrate_fun(x,b,c,R,po,poprime,rhoinfo,pr,s,ps,rhoint);
%         pip(i,j)=r;
%     end
% end
% 
% %%
% sign=[-1 1]; 
% 
% podelta=.001;
% eps=podelta/2;
% maxcount=100;
% po_init=.9;
% postore=zeros(maxcount+1,1);
% postore(1)=po_init;
% ponow=po_init;
% xstore=zeros(maxcount+1,3);
% count=0;
% while count<=maxcount
% 
%     e=equilibria_assortment_interactions(b,c,R,ponow,rhoinfo,pr,s,ps,rhoint);
%     x=e(1,:);
%     x1=x(1);
%     x2=x(2);
%     x3=x(3);
%     x4=0;
%     xstore(count+1,:)=x;
%     
%     pograd=pogradient(x,b,c,R,ponow,rhoinfo,pr,s,ps,rhoint);
%     ponow=ponow+sgn(pograd)*podelta;
%     ponow=round(ponow*1000)/1000;
%     if ponow>=(0-eps) && ponow<=(1+eps)
%         count=count+1;   
%         postore(count+1)=ponow;
%     else
%         postore((count+2):end)=[];
%         xstore((count+2):end,:)=[];
%         count=maxcount+1;
%     end        
% end

