function X = interaction_matrix_round(x,rho)
    x=round(x*100000)/100000;
    X=zeros(3,3);
    X(3,1)=x(1)/(x(1)+x(2)+rho*x(3));
    X(3,2)=x(2)/(x(1)+x(2)+rho*x(3));
    X(3,3)=rho*x(3)/(x(1)+x(2)+rho*x(3));
    X(1,3)=x(3)/(x(1)+x(2)+rho*x(3));
    X(2,3)=x(3)/(x(1)+x(2)+rho*x(3));
    X(1,1)=x(1)/(x(1)+x(2))*(1-X(1,3));
    X(1,2)=x(2)/(x(1)+x(2))*(1-X(1,3));
    X(isnan(X))=0;
%     X(1,:)=X(1,:)/sum(X(1,:));
    X(2,:)=X(1,:);
    
end