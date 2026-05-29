# Bayesian Two-Phase Analysis
# General Model with no interaction terms
# SRS at Phase I

source("tps.q")

# M is dimension of exposure variable X (m_x in paper)
# K is dimension of confounder variable Z (m_z in paper)
# nvec is 2 x m_x x m_z table of phase II data
# Nvec is vector of length 2m_z of phase I sample sizes

# function to perform Bayesian analysis
Bayes2phase.func<-function(Nvec, nvec, maxiter) {
	M<-dim(nvec)[1]
	K<-dim(nvec)[2]/2

	# function to calculate joint P(Y=y,X=x,Z=z) probabilities from lambda
	lambdaprob<-function(lambdas) {
		lambda.yvec<-c(lambdas[1],-lambdas[1])
		lambda.xvec<-c(lambdas[2:M],-sum(lambdas[2:M]))
		lambda.zvec<-c(lambdas[(M+1):(K+M-1)],-sum(lambdas[(M+1):(K+M-1)]))
		lambda.yxmat<-rbind(c(lambdas[(M+K):(2*M+K-2)],-sum(lambdas[(M+K):(2*M+K-2)])),-c(lambdas[(M+K):(2*M+K-2)],-sum(lambdas[(M+K):(2*M+K-2)])))
		lambda.yzmat<-rbind(c(lambdas[(2*(M-1)+K+1):(2*(M-1)+2*(K-1)+1)],-sum(lambdas[(2*(M-1)+K+1):(2*(M-1)+2*(K-1)+1)])),-c(lambdas[(2*(M-1)+K+1):(2*(M-1)+2*(K-1)+1)],-sum(lambdas[(2*(M-1)+K+1):(2*(M-1)+2*(K-1)+1)])))
		lambda.xzmat.temp1<-matrix(lambdas[(2*(M-1)+2*(K-1)+2):(2*(M-1)+2*(K-1)+(M-1)*(K-1)+1)], nrow=M-1, ncol=K-1)
		lambda.xzmat.temp2<-cbind(lambda.xzmat.temp1, -apply(lambda.xzmat.temp1,1,sum))
		lambda.xzmat<-rbind(lambda.xzmat.temp2, -apply(lambda.xzmat.temp2, 2, sum))
	
		qmat.t1<-matrix(-999, nrow=M, ncol=K) # Y=0
		qmat.t2<-matrix(-999, nrow=M, ncol=K) # Y=1
		for(k in 1:M) {
			for(j in 1:K) {
				qmat.t1[k,j]<-exp(lambda.yvec[1]+lambda.xvec[k]+lambda.zvec[j]+lambda.yxmat[1,k]+lambda.yzmat[1,j]+lambda.xzmat[k,j]) #Y=0
				qmat.t2[k,j]<-exp(lambda.yvec[2]+lambda.xvec[k]+lambda.zvec[j]+lambda.yxmat[2,k]+lambda.yzmat[2,j]+lambda.xzmat[k,j]) #Y=1
			}
		}
		qmat1<-qmat.t1/sum(qmat.t1+qmat.t2) # P(Y=0,X=x,Z=z)
		qmat2<-qmat.t2/sum(qmat.t1+qmat.t2) # P(Y=1,X=x,Z=z)
		qmat<-matrix(-999, nrow=M, ncol=2*K)
		for(i in 1:K) {
			qmat[,2*i-1]<-qmat1[,i]
			qmat[,2*i]<-qmat2[,i]
		}
		return(qmat)
	}

	
	# calculate total population size
	N<-sum(Nvec)

	# calculate phase II sample sizes from phase II data
	ntot<-apply(nvec, 2, sum)

	y<-rep(rep(c(0,1),K), ntot)
	x<-rep(rep(c(0:(M-1)),2*K), as.vector(nvec))
	z<-rep(c(0:(K-1)), apply(cbind(ntot[seq(1,2*K,by=2)],ntot[seq(2,2*K,by=2)]),1,sum))

	# fit NPML method to get starting values and proposal variance-covariance matrix
	options(contrasts=c("contr.treatment", "contr.poly"))
	gr<-z+1
	fitML<-tps(y~factor(x)+factor(z), nn0=Nvec[seq(1,2*K,by=2)], nn1=Nvec[seq(2,2*K,by=2)], group=gr, method="ML")

	YXZ.estmat<-matrix(-999, nrow=M, ncol=2*K)
	for(j in 1:M) {
		YXZ.estmat[j,]<-round((nvec[j,]/ntot)*Nvec)
	}
	options(contrasts=c("contr.sum", "contr.poly")) # sum-to-0 constraints
	case.poi<-rep(rep(c(1,2),each=M), K) # 2=case
	xpoi<-rep(c(1:M), 2*K)
	zpoi<-rep(c(1:K),each=2*M)
	countpoi<-as.vector(YXZ.estmat)
	modsat<-glm(countpoi~as.factor(case.poi)*as.factor(xpoi)+as.factor(case.poi)*as.factor(zpoi)+as.factor(xpoi)*as.factor(zpoi), family=poisson)
	varcov.est<-vcov(modsat)[-1,-1]

	# variance-covariance matrix for prior distribution on theta^Y
	# build (m.x-1)+(m.z-1)+1 transformation matrix C
	C.mat<-matrix(0, nrow=(M-1)+(K-1)+1, ncol=(M-1)+(K-1)+1)
	C.mat[1,]<-c(-2,-2,rep(0,M-2),-2,rep(0,K-2))
	C.mat[,2]<-c(-2,rep(2,M-2),4,rep(0,K-1))
	C.mat[M,]<-c(0,4,rep(2,M-2),rep(0,K-1))
	C.mat[,M+1]<-c(-2,rep(0,M-1),rep(2,K-2),4)
	C.mat[M+K-1,]<-c(rep(0,M),4,rep(2,K-2))
	if(M==2 & K!=2) {
		ifelse(is.matrix(C.mat[(M+1):(M+K-2),(M+2):(M+K-1)]),diag(C.mat[(M+1):(M+K-2),(M+2):(M+K-1)])<--2, C.mat[(M+1):(M+K-2),(M+2):(M+K-1)]<--2)
	} else if(M!=2 & K==2) {
		ifelse(is.matrix(C.mat[2:(M-1),3:M]), diag(C.mat[2:(M-1),3:M])<--2, C.mat[2:(M-1),3:M]<--2)
	} else if(M!=2 & K!=2) {
		ifelse(is.matrix(C.mat[2:(M-1),3:M]), diag(C.mat[2:(M-1),3:M])<--2, C.mat[2:(M-1),3:M]<--2)
		ifelse(is.matrix(C.mat[(M+1):(M+K-2),(M+2):(M+K-1)]),diag(C.mat[(M+1):(M+K-2),(M+2):(M+K-1)])<--2, C.mat[(M+1):(M+K-2),(M+2):(M+K-1)]<--2)
	}
	# variance-covariance matrix for betas: uninformative priors
	Sigma.mat<-matrix(0, nrow=(M-1)+(K-1)+1, ncol=(M-1)+(K-1)+1)
	diag(Sigma.mat)<-100^2
	Sigmat<-solve(C.mat)%*%Sigma.mat%*%t(solve(C.mat))

	eyeM<-matrix(0, nrow=M-1, ncol=M-1) # (M-1)x(M-1) identity matrix
	diag(eyeM)<-1
	eyeK<-matrix(0, nrow=K-1, ncol=K-1) # (K-1)x(K-1) identity matrix
	diag(eyeK)<-1
	J.M<-matrix(1, nrow=M-1, ncol=M-1) # (M-1)x(M-1) matrix of 1's
	J.K<-matrix(1, nrow=K-1, ncol=K-1) # (K-1)x(K-1) matrix of 1's

	A<-eyeM-(1/M)*J.M
	B<-eyeK-(1/K)*J.K

	# tuning the constant parameter
	# start with const=0.5
	const<-0.5
	maxcon<-50000
	lambda.mat<-modsat$coef[-1]
	accept.lambdas<-0
	for(i in 1:maxcon) {
		prop.lambdas<-rmvnorm(1, lambda.mat, const*varcov.est)
	
		######################################################
		# accept/reject lambdas using Metropolis-Hastings step
		######################################################
		# P(Y=y,X=x,Z=z)
		qmat<-lambdaprob(lambda.mat) # current lambdas
		qmat.p<-lambdaprob(prop.lambdas) # proposed lambdas
		# P(Y=y,Z=z)
		p.xz.vec<-apply(qmat,2,sum) # current lambdas
		p.xz.vec.p<-apply(qmat.p,2,sum) # proposed lambdas
		
		# P(Ny.z|lambda) for proposed lambdas
		lnumer1.p<-dmultinom(Nvec, N, prob=as.vector(p.xz.vec.p), log=TRUE)
	
		# P(Ny.z|lambda) for lambdas
		ldenom1.p<-dmultinom(Nvec, N, prob=as.vector(p.xz.vec), log=TRUE)
	
		# P(nyxz|lambda) for proposed lambdas
		lnumer1.pn<-rep(-999,2*K)
		for(j in 1:(2*K)) {
			lnumer1.pn[j]<-dmultinom(nvec[,j], ntot[j], prob=qmat.p[,j]/p.xz.vec.p[j], log=TRUE)
		}

		# P(nyxz|lambda) for lambdas
		ldenom1.pn<-rep(-999,2*K)
		for(j in 1:(2*K)) {
			ldenom1.pn[j]<-dmultinom(nvec[,j], ntot[j], prob=qmat[,j]/p.xz.vec[j], log=TRUE)
		}
		
		# normal priors on lambdas
		lnum.priors.p<-rep(-999,4) # 4 groups: (Y,YX,YZ),X,Z,XZ
		ldenom.priors.p<-rep(-999,4) # 4 groups: (Y,YX,YZ),X,Z,XZ	
		# prior for (lambda_y, lambda_yx, lambda_yz)
		lnum.priors.p[1]<-dmvnorm(prop.lambdas[c(1,(M+K):(2*M+K-2),(2*(M-1)+K+1):(2*(M-1)+2*(K-1)+1))], rep(0, K+M-1), Sigmat, log=TRUE) # prior for betas
		ldenom.priors.p[1]<-dmvnorm(lambda.mat[c(1,(M+K):(2*M+K-2),(2*(M-1)+K+1):(2*(M-1)+2*(K-1)+1))], rep(0, K+M-1), Sigmat, log=TRUE)
		# prior for lambda_x 
		lnum.priors.p[2]<-dmvnorm(prop.lambdas[2:M], rep(0, M-1), M*A, log=TRUE)
		ldenom.priors.p[2]<-dmvnorm(lambda.mat[2:M], rep(0, M-1), M*A, log=TRUE)
		# prior for lambda_z 
		lnum.priors.p[3]<-dmvnorm(prop.lambdas[(M+1):(K+M-1)], rep(0, K-1), K*B, log=TRUE)
		ldenom.priors.p[3]<-dmvnorm(lambda.mat[(M+1):(K+M-1)], rep(0, K-1), K*B, log=TRUE)
		# prior for lambda_xz
		lnum.priors.p[4]<-dmvnorm(prop.lambdas[(2*(M-1)+2*(K-1)+2):(2*(M-1)+2*(K-1)+(M-1)*(K-1)+1)], rep(0,(M-1)*(K-1)), (M*(A%x%B))%*%(K*(A%x%B)), log=TRUE)
		ldenom.priors.p[4]<-dmvnorm(lambda.mat[(2*(M-1)+2*(K-1)+2):(2*(M-1)+2*(K-1)+(M-1)*(K-1)+1)], rep(0,(M-1)*(K-1)), (M*(A%x%B))%*%(K*(A%x%B)), log=TRUE)
		
		lnumer.p<-lnumer1.p+sum(lnumer1.pn)+sum(lnum.priors.p)
		ldenom.p<-ldenom1.p+sum(ldenom1.pn)+sum(ldenom.priors.p)
	
		accept_ratio.p<-exp(lnumer.p-ldenom.p)
		if(runif(1)<accept_ratio.p) {
			lambda.mat<-prop.lambdas
			accept.lambdas<-accept.lambdas+1
		} else {
			lambda.mat<-lambda.mat
		}
		if(i%%5000==0) { # every 5000 iterations, we tune const
			acc.rate<-accept.lambdas/5000
			if(acc.rate>0.30) { # if acceptance rate is higher than 30%, increase const
				const<-const+(acc.rate-0.30)*0.80
				#print(c(accept.lambdas, const))
				accept.lambdas<-0
			}
			else { # if acceptance rate is lower than 30%, decrease const
				const<-const-(0.30-acc.rate)*0.80
				#print(c(accept.lambdas, const))
				accept.lambdas<-0
			}
		}
	} # tuning over

	beta.mat<-matrix(-999, nrow=maxiter, ncol=length(fitML$coef), byrow=TRUE)
	lambda.mat<-matrix(rep(modsat$coef[-1],maxiter+1), ncol=dim(varcov.est)[1],byrow=TRUE)
	accept.lambdas<-0

		for(i in 1:maxiter) {
		prop.lambdas<-rmvnorm(1, lambda.mat[i,], const*varcov.est)
	
		######################################################
		# accept/reject lambdas using Metropolis-Hastings step
		######################################################
		# P(Y=y,X=x,Z=z)
		qmat<-lambdaprob(lambda.mat[i,]) # current lambdas
		qmat.p<-lambdaprob(prop.lambdas) # proposed lambdas
		# P(Y=y,Z=z)
		p.xz.vec<-apply(qmat,2,sum) # current lambdas
		p.xz.vec.p<-apply(qmat.p,2,sum) # proposed lambdas
		
		# P(Ny.z|lambda) for proposed lambdas
		lnumer1.p<-dmultinom(Nvec, N, prob=as.vector(p.xz.vec.p), log=TRUE)
	
		# P(Ny.z|lambda) for lambdas
		ldenom1.p<-dmultinom(Nvec, N, prob=as.vector(p.xz.vec), log=TRUE)
	
		# P(nyxz|lambda) for proposed lambdas
		lnumer1.pn<-rep(-999,2*K)
		for(j in 1:(2*K)) {
			lnumer1.pn[j]<-dmultinom(nvec[,j], ntot[j], prob=qmat.p[,j]/p.xz.vec.p[j], log=TRUE)
		}

		# P(nyxz|lambda) for lambdas
		ldenom1.pn<-rep(-999,2*K)
		for(j in 1:(2*K)) {
			ldenom1.pn[j]<-dmultinom(nvec[,j], ntot[j], prob=qmat[,j]/p.xz.vec[j], log=TRUE)
		}
		
		# normal priors on lambdas
		lnum.priors.p<-rep(-999,4) # 4 groups: (Y,YX,YZ),X,Z,XZ
		ldenom.priors.p<-rep(-999,4) # 4 groups: (Y,YX,YZ),X,Z,XZ	
		# prior for (lambda_y,lambda_yx,lambda_yz)
		lnum.priors.p[1]<-dmvnorm(prop.lambdas[c(1,(M+K):(2*M+K-2),(2*(M-1)+K+1):(2*(M-1)+2*(K-1)+1))], rep(0,M+K-1), Sigmat, log=TRUE) # prior for betas
		ldenom.priors.p[1]<-dmvnorm(lambda.mat[i,c(1,(M+K):(2*M+K-2),(2*(M-1)+K+1):(2*(M-1)+2*(K-1)+1))], rep(0,M+K-1), Sigmat, log=TRUE)
		# prior for lambda_x 
		lnum.priors.p[2]<-dmvnorm(prop.lambdas[2:M], rep(0, M-1), M*A, log=TRUE)
		ldenom.priors.p[2]<-dmvnorm(lambda.mat[i,2:M], rep(0, M-1), M*A, log=TRUE)
		# prior for lambda_z 
		lnum.priors.p[3]<-dmvnorm(prop.lambdas[(M+1):(K+M-1)], rep(0, K-1), K*B, log=TRUE)
		ldenom.priors.p[3]<-dmvnorm(lambda.mat[i,(M+1):(K+M-1)], rep(0, K-1), K*B, log=TRUE)
		# prior for lambda_xz
		lnum.priors.p[4]<-dmvnorm(prop.lambdas[(2*(M-1)+2*(K-1)+2):(2*(M-1)+2*(K-1)+(M-1)*(K-1)+1)], rep(0,(M-1)*(K-1)), (M*(A%x%B))%*%(K*(A%x%B)), log=TRUE)
		ldenom.priors.p[4]<-dmvnorm(lambda.mat[i,(2*(M-1)+2*(K-1)+2):(2*(M-1)+2*(K-1)+(M-1)*(K-1)+1)], rep(0,(M-1)*(K-1)), (M*(A%x%B))%*%(K*(A%x%B)), log=TRUE)
		
		lnumer.p<-lnumer1.p+sum(lnumer1.pn)+sum(lnum.priors.p)
		ldenom.p<-ldenom1.p+sum(ldenom1.pn)+sum(ldenom.priors.p)
	
		accept_ratio.p<-exp(lnumer.p-ldenom.p)
		if(runif(1)<accept_ratio.p) {
			lambda.mat[i+1,]<-prop.lambdas
			accept.lambdas<-accept.lambdas+1
		} else {
			lambda.mat[i+1,]<-lambda.mat[i,]
		}
	
		beta0<--2*(lambda.mat[i+1,1]+lambda.mat[i+1,M+K]+lambda.mat[i+1,2*M+K-1])
		if(M==2) {
			betax<-4*lambda.mat[i+1,M+K]
		} else {
			betax<--2*c(lambda.mat[i+1,(M+K+1):(2*M+K-2)]-lambda.mat[i+1,M+K], -sum(lambda.mat[i+1,(M+K):(2*M+K-2)])-lambda.mat[i+1,M+K])
		}		
		if(K==2) {
			betaz<-4*lambda.mat[i+1,2*M+K-1]
		} else {
			betaz<--2*c(lambda.mat[i+1,(2*M+K):(2*(M-1)+2*(K-1)+1)]-lambda.mat[i+1,2*M+K-1], -sum(lambda.mat[i+1,(2*M+K-1):(2*(M-1)+2*(K-1)+1)])-lambda.mat[i+1,2*M+K-1])
		}
		beta.mat[i,]<-c(beta0,betax,betaz)

	} # end of maxiter loop...

	llim<-0.5*maxiter
	betas.m<-apply(beta.mat[c(llim:maxiter),], 2, mean)
	betas.ci<-apply(beta.mat[c(llim:maxiter),], 2, quantile, probs=c(0.025, 0.975))
	return(list(coef=betas.m, CI=t(betas.ci) ))
}
