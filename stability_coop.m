function V = stability_coop(b,c,R,po,rhoinfo,pr,s,ps,rhoint)
    e=equilibria_assortment_interactions(b,c,R,po,rhoinfo,pr,s,ps,rhoint);
    e(2,:)=[];
    e=round(e*10000)/10000;
    num=size(e,1);
    V=zeros(2,num);
    for i=1:num
        J=myjacobian(transpose(e(i,:)),b,c,R,po,rhoinfo,pr,s,ps,rhoint);
        [~,vals]=eig(J);
        V(:,i)=diag(vals);
    end
end