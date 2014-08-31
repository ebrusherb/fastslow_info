function stabeq = stabeq_fun(b,c,R,po,rhoinfo,pr,s,ps,rhoint)
stabeq=zeros(3,1);
e=equilibria_assortment_interactions(b,c,R,po,rhoinfo,pr,s,ps,rhoint);
e(2,:)=[];
e=round(e*10000)/10000;
%     e=[e;1 0 0; 0 1 0; 0 0 1];
V=stability_coop(b,c,R,po,rhoinfo,pr,s,ps,rhoint);
V=real(V);
stable=(sum(V<=0,1)==2);
if sum(stable)>0
    stabeq=e(find(stable>0,1,'first'),:);
end