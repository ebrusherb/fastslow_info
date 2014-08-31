%figuring out why there are gradients even when the costs are the same

%%
povals=0.2:.2:1;
prvals=0.2:.2:1;
postep=.01;
prstep=.01;
Npo=length(povals);
Npr=length(prvals);
% payoff_mat=zeros(Npo,Npr,4);
selection_components=zeros(Npo,Npr,2,4);
coopfreq=zeros(Npo,Npr);
kmat=zeros(Npo,Npr);
gmat=zeros(Npo,Npr);

for i=1:Npo
    for j=1:Npr
        po=povals(i);
        pr=prvals(j);
        
        x=stabeq_fun(b,c,R,po,rhoinfo,pr,s,ps,rhoint);
        coopfreq(i,j)=x(1);
        if sum(x==0)<3
            po1=po/(x(1)+x(2)+rhoinfo*x(3));
            po2=po/(x(1)+x(2)+rhoinfo*x(3));
            po3=po*rhoinfo/(x(1)+x(2)+rhoinfo*x(3));

            k1=kbar(R,po1,pr); 
            k2=kbar(R,po2,pr);
            k3=kbar(R,po3,pr);
            k=k1*x(1)+k2*x(2)+k3*x(3);
            kmat(i,j)=k;

            g=gbar_assortment(R,po,rhoinfo,pr,x(1),x(3),ps,rhoint);
            gmat(i,j)=g;

            grad=poprgradient(x,b,c,R,po,rhoinfo,pr,s,ps,rhoint);
            poprime=po+sgn(grad(1))*postep;
            prprime=pr;
            P=compare_payoffs(x,b,c,R,po,poprime,rhoinfo,pr,prprime,s,ps,rhoint);
            diff=P(:,2)-P(:,1);
            selection_components(i,j,1,:)=sgn(grad(1))*(diff/sum(abs(diff)));
            poprime=po;
            prprime=pr+sgn(grad(2))*prstep;
            P=compare_payoffs(x,b,c,R,po,poprime,rhoinfo,pr,prprime,s,ps,rhoint);
            diff=P(:,2)-P(:,1);
            selection_components(i,j,2,:)=sgn(grad(2))*(diff/sum(abs(diff)));
        end        
%         P = payoffs_divided(x,b,c,R,po,rhoinfo,pr,s,ps,rhoint);
%         payoff_mat(i,j,:)=P;
    end
end
%%
m=min(min(min(min(selection_components))));
M=max(max(max(max(selection_components))));
%%
figsize=[3 4];

figure
subplot(figsize(1),figsize(2),sub2ind(fliplr(figsize),1,1))
imagesc(povals,prvals,transpose(selection_components(:,:,1,1)))
set(gca,'ydir','normal')
caxis manual
caxis([m M]);
title('Costs to cooperators')
colorbar

subplot(figsize(1),figsize(2),sub2ind(fliplr(figsize),2,1))
imagesc(povals,prvals,transpose(selection_components(:,:,1,2)))
set(gca,'ydir','normal')
caxis manual
caxis([m M]);
title('Costs to strangers')
colorbar

subplot(figsize(1),figsize(2),sub2ind(fliplr(figsize),4,1))
imagesc(povals,prvals,transpose(selection_components(:,:,1,4)))
set(gca,'ydir','normal')
caxis manual
caxis([m M]);
title('Costs to discriminators')
colorbar

subplot(figsize(1),figsize(2),sub2ind(fliplr(figsize),3,1))
imagesc(povals,prvals,transpose(selection_components(:,:,1,3)))
set(gca,'ydir','normal')
caxis manual
caxis([m M]);
title('Benefit from being thought to be good')
colorbar

subplot(figsize(1),figsize(2),sub2ind(fliplr(figsize),1,2))
imagesc(povals,prvals,transpose(selection_components(:,:,2,1)))
set(gca,'ydir','normal')
caxis manual
caxis([m M]);
title('Costs to cooperators')
colorbar

subplot(figsize(1),figsize(2),sub2ind(fliplr(figsize),2,2))
imagesc(povals,prvals,transpose(selection_components(:,:,2,2)))
set(gca,'ydir','normal')
caxis manual
caxis([m M]);
title('Costs to strangers')
colorbar

subplot(figsize(1),figsize(2),sub2ind(fliplr(figsize),3,2))
imagesc(povals,prvals,transpose(selection_components(:,:,2,3)))
set(gca,'ydir','normal')
caxis manual
caxis([m M]);
title('Benefit from being thought to be good')
colorbar

subplot(figsize(1),figsize(2),sub2ind(fliplr(figsize),4,2))
imagesc(povals,prvals,transpose(selection_components(:,:,2,4)))
set(gca,'ydir','normal')
caxis manual
caxis([m M]);
title('Costs to discriminators')
colorbar

subplot(figsize(1),figsize(2),sub2ind(fliplr(figsize),1,3))
imagesc(povals,prvals,transpose(coopfreq))
set(gca,'ydir','normal')
title('Cooperators')
colorbar

subplot(figsize(1),figsize(2),sub2ind(fliplr(figsize),2,3))
imagesc(povals,prvals,transpose(kmat))
set(gca,'ydir','normal')
title('K')
colorbar

subplot(figsize(1),figsize(2),sub2ind(fliplr(figsize),3,3))
imagesc(povals,prvals,transpose(gmat))
set(gca,'ydir','normal')
title('G')
colorbar

%%
totalindices=Npo*Npr*2*3;

figure
subplot(2,4,3)
bigindices=sub2ind([Npo Npr 2 3],rows,cols,1*ones(size(rows)),3*ones(size(rows)));
plot(col(coopfreq),col(selection_components(:,:,1,3)),'o')
hold on
plot(col(coopfreq(smallindices)),col(selection_components(bigindices)),'or','MarkerSize',5)
xlabel('Coop freq');ylabel('Benefit from being good');

subplot(2,4,1)
bigindices=sub2ind([Npo Npr 2 3],rows,cols,1*ones(size(rows)),1*ones(size(rows)));
plot(col(kmat),col(selection_components(:,:,1,1)),'o')
hold on
plot(col(kmat(smallindices)),col(selection_components(bigindices)),'or','MarkerSize',5)
xlabel('Knowledge');ylabel('Costs to coops');

subplot(2,4,2)
bigindices=sub2ind([Npo Npr 2 3],rows,cols,1*ones(size(rows)),2*ones(size(rows)));
plot(col(kmat),col(selection_components(:,:,1,2)),'o')
hold on
plot(col(kmat(smallindices)),col(selection_components(bigindices)),'or','MarkerSize',5)
xlabel('Knowledge');ylabel('Costs to strangers');

subplot(2,4,4)
bigindices=sub2ind([Npo Npr 2 3],rows,cols,1*ones(size(rows)),3*ones(size(rows)));
plot(col(gmat),col(selection_components(:,:,1,3)),'o')
hold on
plot(col(gmat(smallindices)),col(selection_components(bigindices)),'or','MarkerSize',5)
xlabel('Good');ylabel('Benefit from being good');

subplot(2,4,7)
bigindices=sub2ind([Npo Npr 2 3],rows,cols,2*ones(size(rows)),3*ones(size(rows)));
plot(col(coopfreq),col(selection_components(:,:,2,3)),'o')
hold on
plot(col(coopfreq(smallindices)),col(selection_components(bigindices)),'or','MarkerSize',5)
xlabel('Coop freq');ylabel('Benefit from being good');

subplot(2,4,5)
bigindices=sub2ind([Npo Npr 2 3],rows,cols,2*ones(size(rows)),1*ones(size(rows)));
plot(col(kmat),col(selection_components(:,:,2,1)),'o')
hold on
plot(col(kmat(smallindices)),col(selection_components(bigindices)),'or','MarkerSize',5)
xlabel('Knowledge');ylabel('Costs to coops');

subplot(2,4,6)
bigindices=sub2ind([Npo Npr 2 3],rows,cols,2*ones(size(rows)),2*ones(size(rows)));
plot(col(kmat),col(selection_components(:,:,2,2)),'o')
hold on
plot(col(kmat(smallindices)),col(selection_components(bigindices)),'or','MarkerSize',5)
xlabel('Knowledge');ylabel('Costs to strangers');

subplot(2,4,8)
bigindices=sub2ind([Npo Npr 2 3],rows,cols,2*ones(size(rows)),3*ones(size(rows)));
plot(col(gmat),col(selection_components(:,:,2,3)),'o')
hold on
plot(col(gmat(smallindices)),col(selection_components(bigindices)),'or','MarkerSize',5)
xlabel('Good');ylabel('Benefit from being good');

%%
figure
plot(col(kmat),col(coopfreq),'o')
hold on
plot(col(kmat(smallindices)),col(coopfreq(smallindices)),'or','MarkerSize',5)
xlabel('Knowledge');ylabel('Coop freq')
%%
v=col(selection_components(:,:,1,2));
v2=col(kmat);
% smallindices=find(v<.19);
% subset=[3 9 10 11];
% smallindices=smallindices(subset);

% smallindices=[21 22 23 24 25 16 17];
% smallindices=[21:25];
% smallindices=16:17;
smallindices=[21 22 23 16 17];
[rows, cols]=ind2sub([Npo Npr],smallindices);
bigindices=sub2ind([Npo Npr 2 3],rows,cols,1*ones(size(rows)),3*ones(size(rows)));

%okay so i can show how selection depends on coop freq, knowledge, and
%knowledge of goodness!  that's great!  and additionally i see that
%something's different when it's an interior equilibrium as opposed to a
%cooperative equilibrium.  good sunday morning work.  