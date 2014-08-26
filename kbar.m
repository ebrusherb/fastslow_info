function K = kbar(R,po,pr)
K=0;
if po~=0
    K=R*pr*po/(1-pr*(1-po))-pr*po*(1-power(pr,R)*power(1-po,R))/power(1-pr*(1-po),2);
end
end