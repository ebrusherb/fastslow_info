library(fields)

#############
############# #functions: 
good.mut<-function(t,po,pr,poprime,prprime,xcoop,xdisc,ps=NULL){
	if(is.null(ps)){ps=.5}
	#g=data.frame("33"=0,"34"=0,"43"=0)
	gold=matrix(0,nrow=3,ncol=1)	
	g=matrix(0,nrow=3,ncol=1)
	if(t>1){
	for(time in 2:t){
		g[1,1]=po*pr*ps+po*pr*(xcoop-ps)*known(time-1,po,pr)+po*pr*xdisc*gold[1,1]+(1-po)*pr*gold[1,1]
		g[2,1]=po*pr*ps+pr*po*(xcoop-ps)*known(time-1,poprime,prprime)+pr*po*xdisc*gold[3,1]+pr*(1-po)*gold[2,1]
		g[3,1]=poprime*prprime*ps+prprime*poprime*(xcoop-ps)*known(time-1,po,pr)+prprime*poprime*xdisc*gold[1,1]+prprime*(1-poprime)*gold[3,1]
		gold=g
		}
		}
	g=matrix(g)
	return(g)
	}

gbar.mut<-function(R,po,pr,poprime,prprime,xcoop,xdisc,ps=NULL){
	t=matrix(1:R,nrow=1)
	g=apply(t,2,good.mut,po=po,pr=pr,poprime=poprime,prprime=prprime,xcoop=xcoop,xdisc=xdisc)
	gbar=apply(g,1,sum)
	return(gbar)
}

payoff_diff<-function(b,c=1,R,po,pr,poprime,prprime,costfun,ps=NULL,disconly=FALSE){
	if(is.null(ps)){ps=.5}
	if(!disconly){
	x=eq2(b,c,R,po,pr,s[1],ps)} else { x=c(0,0,1)}
	diff=NA
	if(x[3]>0){
	k=kbar(R,po,pr)
	kprime=kbar(R,poprime,prprime)
	g=gbar.mut(R,po,pr,poprime,prprime,x[1],x[3],ps)
	g33=g[1]
	g34=g[2]
	g43=g[3]
	p3=-c(R-k)*ps+(b*R-c*k)*x[1]+((b-c)*g33+b*(R-k)*ps)*x[3]-costfun(po,pr)
	p4=-c(R-kprime)*ps+(b*R-c*kprime)*x[1]+(b*g34-c*g43+b*(R-k)*ps)*x[3]-costfun(poprime,prprime)
	diff=p4-p3
	diff=c*(ps-x[1])*(kprime-k)+(b*g34-c*g43-(b-c)*g33)*x[3]-costfun(poprime,prprime)+costfun(po,pr)
	}
	return(diff)
}

cost<-function(po,pr){
	#synergy=0
	#pocost=.0001
	#prcost=.0001
	# pocost=0
	# prcost=0
	constant=0.01
	tot=synergy*po*pr+pocost*po+prcost*pr+constant
	#tot=.01
	return(tot)
}

selection<-function(b,c=1,R,po,pr,costfun,ps=NULL,disconly=FALSE,eps=.001){
	if(is.null(ps)){ps=.5}
		if(po==1){
			po_down=payoff_diff(b,c,R,po,pr,po-eps,pr,costfun,ps,disconly)
			dpo=-po_down/eps
		} else if(po==0){
			po_up=payoff_diff(b,c,R,po,pr,po+eps,pr,costfun,ps,disconly)
			dpo=po_up/eps
		} else{
		po_up=payoff_diff(b,c,R,po,pr,po+.5*eps,pr,costfun,ps,disconly)
		po_down=payoff_diff(b,c,R,po,pr,po-.5*eps,pr,costfun,ps,disconly)
		dpo=(po_up-po_down)/(eps)
		} 
		if(pr==1){
			pr_down=payoff_diff(b,c,R,po,pr,po,pr-eps,costfun,ps,disconly)
			dpr=-pr_down/eps
		} else if(pr==0){
			pr_up=payoff_diff(b,c,R,po,pr,po,pr+eps,costfun,ps,disconly)
			dpr=pr_up/eps
		} else{
		pr_up=payoff_diff(b,c,R,po,pr,po,pr+.5*eps,costfun,ps,disconly)
		pr_down=payoff_diff(b,c,R,po,pr,po,pr-.5*eps,costfun,ps,disconly)
		dpr=(pr_up-pr_down)/(eps)
		} 
	return(c(dpo,dpr))
}

# ##########################
# ####################### # creat pairwise invasibility plot
# #### selection on po
# R=10
# b=10
# c=1
# s=0.01
# ps=.5

# l=10
# povals=seq(.01,1,length.out=l)

# pr=.5

# diffvals=array(1,c(l,l))

# for(i in 1:l){
	# for(j in 1:l){
	# pores=povals[i]
	# pomut=povals[j]
	# d=payoff_diff(b,c,R,pores,pr,pomut,pr,c(s,s))
	# diffvals[i,j]=d
	# }
# }

# res=povals[which(!is.na(apply(diffvals,1,max)))]
# diffvals=diffvals[which(!is.na(apply(diffvals,1,max))),]

# # filename=paste("PIP_" , "b=" , round(b,2), "_R=", R, "_pr=", pr, "_s=", s, ".png",sep="")
# # png(file=filename)
# M=max(abs(range(diffvals)))
# v=seq(-M,M,length.out=21)
# xlab=bquote("Resident" ~ p[o])
# ylab=bquote("Mutatn" ~ p[o])
# image(res,povals,diffvals,breaks=v,col=colorRampPalette(c("red","white","blue"))(length(v)-1),xlab=xlab,ylab=ylab)
# image.plot(res,povals,diffvals,breaks=v,col=colorRampPalette(c("red","white","blue"))(length(v)-1),xlab=xlab,ylab=ylab,axes=FALSE)
# box()
# abline(0,1,lwd=1)
# title("Mutant Growth Rate")
# mtext(bquote(b * "=" * .(round(b,2)) * "," ~ "R = " * .(R) * "," ~ p[r] * "=" * .(pr) *  "," ~ s * "=" * .(s)))

# L=contourLines(res,povals,diffvals,levels=0)
# m=length(L)

# est=NULL
# isocline=NULL
# for(i in 1:m){
	# c=L[[i]]
	# slope=diff(c$y)/diff(c$x)
	# find=data.frame(slope=slope,diff=abs(c$y-c$x)[-1],pox=c$x[-1],poy=c$y[-1])
	# onthediag=which(abs(slope-1)<.5)
	# onthediag=which(slope>.1)
	# if(length(onthediag)>0){
		# find=find[-onthediag,]}
	# est=c(est,find[,3][which.min(find[,2])])
	# isocline=rbind(isocline,data.frame(x=find$pox,y=find$poy))
# }

# est=mean(est)
# isocline=isocline[order(isocline[,1]),]

# points(est,est,cex=2,lwd=2)

# lines(isocline$x,isocline$y,lwd=1)
# contour(res,povals,diffvals,levels=0,add=TRUE)
# axis(1,at=sort(c(round(est,2),seq(.2,1,by=.2))),labels=sort(c(round(est,2),seq(.2,1,by=.2))))
# axis(2,at=sort(c(round(est,2),seq(.2,1,by=.2))),labels=sort(c(round(est,2),seq(.2,1,by=.2))))

# # dev.off()

# ########
# ##selection on pr


# R=10
# b=5
# c=1
# s=0.01
# ps=.75

# synergy=0
# pocost=.5
# prcost=.5

# #povals=seq(0.1,1,length.out=9)
# #po=.5
# povals=c(0,.1,.5,1)

# #windows()
# layout(matrix(1:4,nrow=2,byrow=TRUE))

# for(o in 1:length(povals)){
# po=povals[o]
# l=50
# prvals=seq(.01,1,length.out=l)

# diffvals=array(1,c(l,l))

# for(i in 1:l){
	# for(j in 1:l){
	# prres=prvals[i]
	# prmut=prvals[j]
	# d=payoff_diff(b,c,R,po,prres,po,prmut,cost,ps,disconly=TRUE)
	# diffvals[i,j]=d
	# }
# }

# res=prvals[which(!is.na(apply(diffvals,1,max)))]
# diffvals=diffvals[which(!is.na(apply(diffvals,1,max))),]

# #windows()
# # filename=paste("PIP_" , "b=" , round(b,2), "_R=", R, "_pr=", pr, "_s=", s, ".png",sep="")
# # png(file=filename)
# M=max(abs(range(diffvals)))
# v=seq(-M,M,length.out=21)
# xlab=bquote("Resident" ~ p[r])
# ylab=bquote("Mutatn" ~ p[r])

# #par(oma=c(0,0,0,1))
# #image(res,prvals,diffvals,breaks=v,col=colorRampPalette(c("red","white","blue"))(length(v)-1),xlab=xlab,ylab=ylab)
# image.plot(res,prvals,diffvals,breaks=v,col=colorRampPalette(c("red","white","blue"))(length(v)-1),xlab=xlab,ylab=ylab,axes=FALSE)
# box()
# abline(0,1)
# title("Mutant Growth Rate")
# #mtext(bquote(b * "=" * .(round(b,2)) * "," ~ "R = " * .(R) * "," ~ p[o] * "=" * .(po) *  "," ~ s * "=" * .(s)))
# mtext(bquote(p[o] * "=" * .(po)))

# L=contourLines(res,prvals,diffvals,levels=0)
# m=length(L)

# est=NULL
# isocline=NULL
# for(i in 1:m){
	# cont=L[[i]]
	# slope=diff(cont$y)/diff(cont$x)
	# find=data.frame(slope=slope,diff=abs(cont$y-cont$x)[-1],pox=cont$x[-1],poy=cont$y[-1])
	# onthediag=which(abs(slope-1)<1)
	# #onthediag=which(slope>.1)
	# if(length(onthediag)>0){
		# find=find[-onthediag,]}
	# est=c(est,find[,3][which.min(find[,2])])
	# isocline=rbind(isocline,data.frame(x=find$pox,y=find$poy))
# }

# est=mean(est)
# isocline=isocline[order(isocline[,1]),]

# points(est,est,cex=2,lwd=2)

# lines(isocline$x,isocline$y,lwd=1)
# #contour(res,prvals,diffvals,levels=0,add=TRUE)
# axis(1,at=sort(c(round(est,2),seq(.2,1,by=.2))),labels=sort(c(round(est,2),seq(.2,1,by=.2))))
# axis(2,at=sort(c(round(est,2),seq(.2,1,by=.2))),labels=sort(c(round(est,2),seq(.2,1,by=.2))))

# # dev.off()
# }

