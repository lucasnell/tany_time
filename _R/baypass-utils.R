###############################################*
###function: Simulate under the inference model
##############################################*
simulate.baypass <- function(omega.mat,
                             nsnp=1000,
                             beta.coef=NA,
                             beta.pi=c(1,1),
                             pop.trait=0,
                             sample.size=100,
                             pi.maf=0.05,
                             suffix="sim",
                             print.sim.params.values=FALSE,
                             remove.fixed.loci=FALSE,
                             output.bayenv.format=FALSE,
                             coverage=NA){
  #sample.size=matrix with npop colums or vector of length npop with count. If matrix and poolseq data, the population sample size (i.e., haploid size of the pools for Pool-Seq data) is set to the colMax of sample.size
  #coverage = matrix with npop colums or vector of length npop with coverage => activate poolseq data
  #pop.trait = vector of length npop OR matrix of dimension npop*ntrait
  #beta.coef = vector or length nsnp_asso OR matrix of dimension nsnp_asso*ntrait

  if(!is.matrix(omega.mat)){stop("omega.mat must be a matrix (e.g., as obtained from *mat_omega.out BayPass output)")}
  npop=nrow(omega.mat)

  if(sum(is.na(beta.coef))==0){
    beta.coef=as.matrix(beta.coef) #security
    pop.trait=as.matrix(pop.trait)
    if(ncol(beta.coef)!=ncol(pop.trait)){
      stop(paste0(ncol(beta.coef)," traits specified in beta.coef while ",ncol(pop.trait)," specified in pop.trait: these must be the same"))
    }
    if(nrow(pop.trait)!=npop){
      stop(paste0(nrow(pop.trait)," specified in pop.trait while ",npop," in omega.mat: these must be the same"))
    }
    simu.cov=TRUE
    ntraits=ncol(beta.coef)
    nsnp.asso=nrow(beta.coef)
    alpha.cov = beta.coef %*% t(pop.trait)
    if(nsnp>0){
      beta.coef=rbind(matrix(0,nsnp,ntraits),beta.coef)
      alpha.cov=rbind(matrix(0,nsnp,npop),alpha.cov)
    }
    nsnp=nsnp + nsnp.asso
  }else{simu.cov=FALSE}

  if(length(sample.size)==1){
    NN=matrix(sample.size,nsnp,npop)
    poolsize=rep(sample.size,npop)
  }else{
    sample.size=as.matrix(sample.size)
    if(ncol(sample.size)==1){#c'est un vecteur
      if(nrow(sample.size)!=npop){stop("Sample size dimension must be of length 1 or have the same size as the rank of omega.mat or must be a matrix with npop columns")}
      NN=matrix(rep(sample.size,nsnp),nsnp,npop,byrow=TRUE)
      poolsize=as.numeric(sample.size)
    }else{
      tmp.snp=sample(1:nrow(sample.size),nsnp,replace=TRUE)
      NN=sample.size[tmp.snp,]
      poolsize=apply(sample.size,2,max)
    }}

  if(length(coverage)==1){
    if(is.na(coverage)){
      poolseq=FALSE
    }else{poolseq=TRUE ; NN.coverage=matrix(coverage,nsnp,npop)}
  }else{
    poolseq=TRUE
    coverage=as.matrix(coverage)
    if(ncol(coverage)==1){#c'est un vecteur
      if(nrow(coverage)!=npop){
        stop("Coverage dimension must be of length 1, or have the same size as the rank of omega.mat or must be a matrix with npop columns")
      }
      NN.coverage=matrix(rep(coverage,nsnp),nsnp,npop,byrow=TRUE)
    }else{
      tmp.snp=sample(1:nrow(coverage),nsnp,replace=TRUE)
      NN.coverage=coverage[tmp.snp,]
    }
  }

  if(poolseq){
    YY.coverage=NN.coverage*0
    NN=matrix(rep(poolsize,nsnp),nsnp,npop,byrow=TRUE)
  }

  if(length(beta.pi)!=2){stop("beta.pi must be of length 2")}else{
    if(sum(beta.pi==1)==2){Pi=runif(nsnp,pi.maf,1-pi.maf)}else{
      Pi=rbeta(nsnp,beta.pi[1],beta.pi[2])
      Pi[Pi<pi.maf]=pi.maf ; Pi[Pi>1-pi.maf]=1-pi.maf
    }
  }

  #pop. allele frequencies
  ALPHA=matrix(rnorm(nsnp*npop),npop,nsnp) #alpha_tilde (transpose to facilitate further rescaling by column)
  ALPHA=t(chol(omega.mat))%*%ALPHA         #transformation into alpha (! lower triangular choleski matrix but chol returns upper triangular)
  ALPHA=t(ALPHA)*sqrt(Pi*(1-Pi)) + Pi
  if(simu.cov){ALPHA=ALPHA + alpha.cov}
  ALPHA_tr=ALPHA ; ALPHA_tr[ALPHA>1]=1 ; ALPHA_tr[ALPHA<0]=0
  #sample allele count
  YY=matrix(0,nsnp,npop)
  for(i in 1:nsnp){
    for(j in 1:npop){
      YY[i,j]=rbinom(1,size=NN[i,j],prob=ALPHA_tr[i,j])
      if(poolseq){
        YY.coverage[i,j]=rbinom(1,size=NN.coverage[i,j],prob=YY[i,j]/NN[i,j])
      }
    }
    if(i%%(nsnp/10)==0){cat(i,"SNP simulated out of",nsnp,"\n")}
  }

  if(output.bayenv.format){remove.fixed.loci=TRUE}#Bayenv doesn't accept fixed loci

  if(remove.fixed.loci){
    if(poolseq){tmp.freq=rowSums(YY.coverage)/rowSums(NN.coverage)}else{tmp.freq=rowSums(YY)/rowSums(NN)}
    snp.sel=tmp.freq>0 & tmp.freq<1
    nsnp=sum(snp.sel)
    YY=YY[snp.sel,] ; NN=NN[snp.sel,] ; Pi=Pi[snp.sel] ; ALPHA=ALPHA[snp.sel,]
    if(poolseq){YY.coverage=YY.coverage[snp.sel,];NN.coverage=NN.coverage[snp.sel,]}
    if(simu.cov){beta.coef=beta.coef[snp.sel]}
    cat("Number of SNPs removed: ",sum(!snp.sel),"\n")
  }

  if(print.sim.params.values){
    write.table(file=paste("pi.",suffix,sep=""),Pi,quote=F,col.names=F,row.names=F)
    write.table(file=paste("alpha.",suffix,sep=""),ALPHA,quote=F,col.names=F,row.names=F)

    if(simu.cov){
      write.table(file=paste("betacoef.",suffix,sep=""),beta.coef,quote=F,col.names=F,row.names=F)
      write.table(file=paste("pheno.",suffix,sep=""),t(pop.trait),quote=F,col.names=F,row.names=F)
    }
  }

  mat_baypass=cbind(YY,NN)
  all2.pos=2*(1:npop)
  mat_baypass[,all2.pos-1]=YY ; mat_baypass[,all2.pos]=NN-YY

  if(!poolseq | print.sim.params.values){
    write.table(file=paste("G.",suffix,sep="") ,mat_baypass,quote=F,col.names=F,row.names=F)
    if(output.bayenv.format){
      mat_bayenv=rbind(YY,NN)
      all2.pos=2*(1:nsnp)
      mat_bayenv[all2.pos-1,]=YY ; mat_bayenv[all2.pos,]=NN-YY
      mat_bayenv=cbind(mat_bayenv,rep("",nsnp))
      write.table(file=paste("bayenv_freq.",suffix,sep=""),mat_bayenv,sep="\t",quote=F,col.names=F,row.names=F)
    }
  }

  if(poolseq){
    mat_baypass=cbind(YY.coverage,NN.coverage)
    all2.pos=2*(1:npop)
    mat_baypass[,all2.pos-1]=YY.coverage ; mat_baypass[,all2.pos]=NN.coverage-YY.coverage
    write.table(file=paste("Gpool.",suffix,sep="") ,mat_baypass,quote=F,col.names=F,row.names=F)

    if(output.bayenv.format){
      mat_bayenv=rbind(YY.coverage,NN.coverage)
      all2.pos=2*(1:nsnp)
      mat_bayenv[all2.pos-1,]=YY.coverage ; mat_bayenv[all2.pos,]=NN.coverage-YY.coverage
      mat_bayenv=cbind(mat_bayenv,rep("",nsnp))
      write.table(file=paste("bayenv_freq_pool.",suffix,sep=""),mat_bayenv,sep="\t",quote=F,col.names=F,row.names=F)
    }
    write.table(file=paste("poolsize.",suffix,sep=""),t(poolsize),quote=F,col.names=F,row.names=F)
  }

  if(simu.cov){
    if(poolseq){
      list(Y.pool=YY.coverage,N.pool=NN.coverage,Y.sim=YY,N.sim=NN,Pi.sim=Pi,alpha.sim=ALPHA,omega.sim=omega.mat,betacoef.sim=beta.coef)
    }else{
      list(Y.sim=YY,N.sim=NN,Pi.sim=Pi,alpha.sim=ALPHA,omega.sim=omega.mat,betacoef.sim=beta.coef)
    }
  }else{
    if(poolseq){
      list(Y.pool=YY.coverage,N.pool=NN.coverage,Y.sim=YY,N.sim=NN,Pi.sim=Pi,alpha.sim=ALPHA,omega.sim=omega.mat)
    }else{
      list(Y.sim=YY,N.sim=NN,Pi.sim=Pi,alpha.sim=ALPHA,omega.sim=omega.mat)
    }
  }
}

###############################################*
###function: transform geno input files into YY and NN
##############################################*

geno2YN<-function(genofile){
 data=read.table(genofile)
 npop=ncol(data)/2 ; all1=seq(1,2*npop,2) ; all2=all1+1
 YY=data[,all1]
 NN=YY+data[,all2]
 list(YY=YY,NN=NN)
}

###############################################*
###function: Compute FMD distance (Forstner and Moonen, 2003) between two covariance matrices
##############################################*
fmd.dist<-function(mat1,
                   mat2){
 require(geigen)
 return(sqrt(sum(log(geigen(mat1,mat2)$values)**2)))
}

###############################################*
###function: Simulate covariate value correlated to a given Omega PC
##############################################*
simulate.PCcorrelated.covariate <- function(omega,
                                            axis=1,
                                            targeted.rho=0.1){
 npops=nrow(omega)
 om.svd=svd(omega)
 PC=om.svd$u[,axis]
 u=scale(PC)/sqrt(npops-1)
 #sample left null space of u by taking the cplt of the projection of a randomly sampled vector
 v.rnd=scale(rnorm(npops))/sqrt(npops-1)
 v.rnd.proj=as.numeric(t(v.rnd)%*%u)*u
 v=v.rnd-v.rnd.proj
 v=v/sqrt(sum(v*v))
 #compute w
 w=targeted.rho*u+sqrt(1-targeted.rho**2)*v
 return(scale(w)[,1])
}

###############################################*
###function: Plot Omega after spectral decomposition: Omega=UDU' with D=diagonal matrix of singular values
##############################################*
plot.omega <- function(omega,
                       PC=c(1,2),
                       pop.names=paste0("Pop",1:nrow(omega)),
                       main=expression("SVD of "*Omega),
                       col=rainbow(nrow(omega)),
                       pch=16,
                       pos=2){
 om.svd=svd(omega)
 eig=om.svd$d
 pcent.var=100*eig/sum(eig)
 plot(om.svd$u[,PC],main=main,pch=pch,col=col,
      xlab=paste0("PC",PC[1]," (",round(pcent.var[PC[1]],2),"%)"),
      ylab=paste0("PC",PC[2]," (",round(pcent.var[PC[2]],2),"%)")
      )
 text(om.svd$u[,PC[1]],om.svd$u[,PC[2]],pop.names,col=col,pos=pos)
 list(PC=om.svd$u,eig=eig,pcent.var=pcent.var)
}

#################################################*
#####function: compute genetic offset from estimated regression coefficient
#################################################*
compute_genetic_offset<-function(beta.coef=NULL,
                                 regfile="summary_betai.out",
                                 candidate.snp=NULL,
                                 covfile="cov.baypass",
                                 newenv=NULL,
                                 refenv=NULL,
                                 scalecov=TRUE,
                                 compute.rona=FALSE){
  #adapted from genetic.gap function of the LEA package to compute genetic gap (Gain and Francois, 2023)
  #beta.coef=matrix of nsnp X ncovariates coefficient. If NULL, user may provide a BayPass output file
  #regfile = BayPass output file (either [outprefix_]summary_betai_reg.out or [outprefix_]summary_betai.out if IS or MC estimates respectively)
  #covfile = BayPass covariate file (-efile)
  #newenv = covariate values for the new environment(s): either a vector of ncovariates is only a single new environment or a matrix of with covariate values for several new environments in column (i.e., ncovariates X nenvironments)
  #refenv = covariable values for the ref environment(s) to be used instead of covfile (i.e.: default: ref environment = those in cov files specified by each row of covariable value). Same format as newenv
  require(data.table)
  ##reading and checking covariate
  cov=as.matrix(read.table(covfile))
  n.cov=nrow(cov) ; n.pop=ncol(cov)
  ##reading and checking newenv
  if(is.null(newenv)){
    stop("Please Provide a vector or a matrix (n_covariates X n_newenv) of covariate values for the target environment(s)\n")
  }else{newenv=as.matrix(newenv)}
  if(nrow(newenv)!=n.cov){
    stop("Check newenv argument: The number of covariates for the new environment is not equal to the number of covariates in the original covariate file\n")
  }
  n.newenv=ncol(newenv)
  ##Scaling pop and new covariates
  if(scalecov){
    m.cov=rowMeans(cov) ; sd.cov=apply(cov,1,sd)
    cov=(cov-m.cov)/sd.cov
    newenv=(newenv-m.cov)/sd.cov
  }
  ##reading and checking refenv
  if(is.null(refenv)){##refenv is set to cov if not provided (default)
    cat("Genetic Offset statistics will be estimated between each of the",n.pop,"population environment (specified by",n.cov,"covariables)\n provided in the original covariate file and the",n.newenv,"target environments provided with newenv argument.\n
  The resulting GO matrix will contain",n.pop,"reference environment (row) times",n.newenv,"target environments (column) entries\n")
    refenv=cov ; n.refenv=n.pop
  }else{
    refenv=as.matrix(refenv)
    if(nrow(refenv)!=n.cov){
      stop("Check refenv argument: the number of covariates for the reference environment is not equal to the number of covariates in the original covariate file\n")
    }
    n.refenv=ncol(refenv)
    ##Scaling pop and new covariates
    if(scalecov){refenv=(refenv-m.cov)/sd.cov}
    cat("Genetic Offset statistics will be estimated between each of the",n.refenv,"reference environment (specified by",n.cov,"covariables)\n provided with refenv argument and the",n.newenv,"target environments provided with newenv argument.\n
  The resulting GO matrix will contain",n.refenv,"reference environment (row) times",n.newenv,"target environments (column) entries\n")
  }
  ##reading regression coefficient
  if(is.null(beta.coef)){
    type=NULL
    if(grepl("summary_betai_reg.out",regfile)){
      type="IS"
      beta.coef=matrix(fread(regfile,data.table=FALSE,header=T)$Beta_is,ncol=n.cov) #nsnp,ncov
    }
    if(grepl("summary_betai.out",regfile)){
      type="MC"
      beta.coef=matrix(fread(regfile,data.table=FALSE,header=T)$M_Beta,ncol=n.cov) #nsnp,ncov
    }
    if(is.null(type)){
      stop("Check regfile argument: No proper BayPass output files with estimates of regression coefficient (i.e., either (either [outprefix_]summary_betai_reg.out or [outprefix_]summary_betai.out if IS or MC estimates respectively) could be found.\nPlease provide a valid BayPass output file (regfile argument) OR a matrix of estimated regression coefficient (beta.coef argument)")
    }
  }else{
    beta.coef=as.matrix(beta.coef)
    if(ncol(beta.coef)!=n.cov){stop("ERROR: The number of covariates (columns) of the beta.coef matrix is not equal to the number of covariates in the original covariate file\n")}
  }

  ###################*
  ggap=matrix(0,n.refenv,n.newenv)
  if(!is.null(candidate.snp)){beta.coef=beta.coef[candidate.snp,]}
  BtB=t(beta.coef)%*%beta.coef/nrow(beta.coef) # cov(beta.coef)
  eig <- eigen(BtB,symmetric = TRUE) #eigendecompositon of BtB=(beta.coef*t(beta.coef))/nsnp=UDU' car symetrique

  for(i in 1:n.refenv){#usually more target than ref
    #    M=t(cov-newenv[,i])%*%t(beta.coef) #sum of reg. coeff for each covariable weighted by the difference between the current and new env. <=> sum. of expected difference in predicted and current allele frequency!
    # E(f_cur_i - f_new_i) = (Pi_i + S_j(beta_ij*Ycur_j)) - (Pi_i + S_j(beta_ij*Ynew_j)) = S_j(beta_ij*(Ycur_j-Ynew_j))
    #    ggap[,i] = rowMeans(M**2) #diag(M %*% t(M))/nrow(beta.coef)
    #    mean.diff.freq[,i]=rowMeans(M)
    diff=newenv-refenv[,i]
    ggap[i,] = colSums(eig$values*(t(eig$vectors)%*%diff)**2)  #GG=(x-x*)%*%BtB%*%(x-x*)'=((x-x*)%*%U)%*%D%*%((x-x*)%*%U)'
  }
  #covariable importance
  covimp=(eig$vectors**2)%*%(eig$values)
  out=list(go = ggap,BtB.eigenvalues = eig$values, BtB.eigenvectors = eig$vectors,covimp=as.vector(covimp))

  if(compute.rona){
    out$rona=matrix(0,n.refenv,n.newenv)
    for(i in 1:n.refenv){out$rona[i,]=rowMeans(abs(t(newenv-refenv[,i])%*%t(beta.coef)))}
  }

  return(out)
}

###############################################*
###function: concatenate sub-datasets results
##############################################*
concatenate_res<-function(dir="./",
                          anaprefix="ana",
                          anasep="_",
                          extension="",
                          nsubsets=2,
                          snpdet_prefix="./detsnp.sub",
                          retrieve_pi_xtx=TRUE,
                          retrieve_bfis=TRUE,
                          retrieve_c2=FALSE){
  #extension should be the same for snpdet files and baypass output files (i.e., all files should be compressed the same way or not)
  #snp_det_prefix should include the path (since it may be different than the directory containing BayPass output files)
  require(data.table)
  corepref=""
  if(nchar(anaprefix)>0){corepref=paste0(anaprefix, anasep)}
  if(nsubsets<1){stop("Check nsubsets: at least one subset is needed\n")}

  for(i in 1:nsubsets){
    # cat("Processing run",i,"out of",nsubsets,"\n")
    tmp.snpdet=fread((paste0(snpdet_prefix,i,extension)),data.table = F)[,1:2]
    colnames(tmp.snpdet)=c("CHR","POS")
    tmp.nsnps=nrow(tmp.snpdet)
    if(retrieve_pi_xtx){
      tmp=fread(paste0(dir,"/",corepref,i,"_summary_pi_xtx.out",extension),data.table=F)[,c("M_P","M_XtX","XtXst")]
      tmp.snpdet=cbind(tmp.snpdet,tmp)
    }
    if(retrieve_bfis){
      bfis_fn <- paste0(dir,"/",corepref,i,"_summary_betai_reg.out",extension)
      if (!file.exists(bfis_fn)) {
          bfis_fn <- paste0(dir,"/",corepref,i,"_summary_betai.out",extension)
          if (!file.exists(bfis_fn)) stop("can't find ", bfis_fn, "!")
      }
      tmp=fread(bfis_fn,data.table=F)$"BF(dB)"
      tmp=matrix(as.numeric(tmp),tmp.nsnps)
      colnames(tmp)=paste0("BFis_cov_",1:ncol(tmp))
      tmp.snpdet=cbind(tmp.snpdet,tmp)
    }
    if(retrieve_c2){
      tmp=fread(paste0(dir,"/",corepref,i,"_summary_contrast.out",extension),data.table=F)$"C2"
      tmp=matrix(tmp,tmp.nsnps)
      colnames(tmp)=paste0("C2_",1:ncol(tmp))
      tmp.snpdet=cbind(tmp.snpdet,tmp)
    }
    if(i==1){res=tmp.snpdet}else{res=rbind(res,tmp.snpdet)}
  }
  res=res[order(res$CHR,res$POS),]
  return(res)
}

####################################*
## function: Compute Local Score and Identify significant windows of markers
####################################*
compute.local.scores<-function(snp.position,                #data.frame with 2 columns: i) chromosome name; ii) position (should be ordered)
                               snp.pi,                      #overall (average) reference allele frequency of the SNPs (=Pi)
                               snp.pvalue=NULL,             #p-value associated to the test statistics (C2 or XtX*) assumed to be on the -log10 scale
                               snp.bf=NULL,                 #BF (in deciban units i.e. 10log10(BF))
                               xi=1,                        #parameter of the lindley process (if not in c(1,2,3,4) local scores pvalues are not computed)
                               pval.local.score.thres=0.01, #local-scores significance thresholds are computed assuming uniform p-values associated to the scores
                               lindley.thr=NULL,            #User defined threshold on lindley score (used to identify windows in place of p-values if xi=1 or xi=2)
                               min.maf=0.2,                 #MAF threshold on Pi (Position with Pi<min.maf are excluded from the computation of the local score)
                               min.nsnp=1e4,                #Minimum number of SNPs in a chromosome (chromosome/scaffolds with <min.nsnp are excluded from the analysis)
                               bf.scale=1,                  #Rescaling factor of the BF (e.g. setting bf.scale to 0.5 halves the BF values in dB)
                               plot.pvalhist=FALSE,         #If TRUE plot the distribution of SNP associated p-values derived from the input statistics
                               manplot=FALSE,               #if TRUE a manhattan plot of the local scores is displayed
                               lindley.capping=NULL,        #Capping value of lindley scores in the manhattan plot
                               ...                          #arguments of the plot function (used when manplot is TRUE)
                               ){
  #Implements computation of local scores as described in Fariello et al. (2017, https://doi.org/10.1111/mec.14141)
  #The code adapted original functions available here: https://forge-dga.jouy.inra.fr/projects/local-score
  #Only values of xi=1, 2, 3 and 4 are considered since these allow for exact computation of confidence thresholds
  #  -for xi=1 or xi=2: see Fariello et al. 2017 for xi=1 and xi=2
  #  -for xi=3 or xi=4: see Bonhomme et al. 2019 (https://www.nature.com/articles/s41437-019-0235-x ; SM)
  #BF are assumed given on dB units. Scores associated to negative BF (i.e., evidence against association) are derived from a uniform draw (i.e., p-value distribution under H0)
  require(data.table)
  #checks
  nsnp=nrow(snp.position)
  chr.names=table(snp.position[,1])
  cat(nsnp," SNPs found in ",length(chr.names)," different chromosomes\n")

  if(length(snp.pi)!=nsnp){stop("snp.pi vector must be the same length as the number of SNPs\n")}
  snp.maf=0.5-abs(0.5-snp.pi)

  if(!is.null(snp.bf)){
    if(!is.null(snp.pvalue)){stop("Please provide either of vectors of BF or a vectors of p-value")}
    if(length(snp.bf)!=nsnp){stop("snp.bf vector must be the same length as the number of SNPs\n")}
    na.remove=is.na(snp.bf)
    if(sum(na.remove)>0){cat("WARNING:",sum(na.remove)," values are NA: the corresponding SNPs will be discarded\n")}
    tmp.inf=is.infinite(snp.pvalue)
    if(sum(tmp.inf)>0){
      cat("WARNING:",sum(tmp.inf)," values are Inf in the BF vector: there were set to the maximum (+Inf) or minimum (-Inf) finite value\n")
      tmp.range=range(snp.bf[is.finite(snp.bf)])
      snp.bf[snp.bf == -Inf & !(na.remove)]=tmp.range[1]
      snp.bf[snp.bf == Inf & !(na.remove)]=tmp.range[2]
    }
    is.bf=TRUE
  }else{
    if(length(snp.pvalue)!=nsnp){stop("snp.pvalue vector must be the same length as the number of SNPs\n")}
    if(min(snp.pvalue,na.rm=T)<0){stop("Some value are negative in snp.pvalue vector: Pvalues must be on -log10 scale\n")}
    if(max(snp.pvalue,na.rm=T)<1){stop("No value >1 in snp.pvalue vector: Pvalues must be on -log10 scale\n")}
    na.remove=is.na(snp.pvalue)
    if(sum(na.remove)>0){cat("WARNING:",sum(na.remove)," values are NA: the corresponding SNPs will be discarded\n")}
    tmp.inf=is.infinite(snp.pvalue)
    if(sum(tmp.inf)>0){
      cat("WARNING:",sum(tmp.inf)," values are Inf in the pvalue vector: there were set to the maximum finite value\n")
      snp.pvalue[tmp.inf]=max(snp.pvalue[is.finite(snp.pvalue)])
    }
    is.bf=FALSE
  }

  for(i in names(chr.names)){
    tmp=snp.position[snp.position[,1]==i,2]
    if(is.unsorted(tmp) & is.unsorted(rev(tmp))){
      stop("WARNING chromosome ",i," SNPs are not sorted (in increasing or decreasing order)\n")
      #  cat("\tThese will be reordered (if chromosome retained)\n")
    }
  }

  snp.sel=snp.maf>min.maf & snp.position[,1]%in%(names(chr.names[chr.names>min.nsnp])) & !(na.remove)
  cat(sum(snp.sel)," SNPs retained mapping to",length(unique(snp.position[snp.sel,1]))," different chromosomes\n")

  if(!is.null(lindley.capping)){
    if(!manplot){
      cat("WARNING: the argument lindley.capping is only used for plotting (i.e., disregarded when manplot is FALSE\n")
    }else{
      if(lindley.capping<0){stop("lindley.capping must be a positive number\n")}
    }
  }

  if(sum(xi != 1:4)==4){
    cat("WARNING: xi is not equal to 1, 2, 3 or 4 => p-values associated with the local-score will not be computed (see Fariello et al., for simulation based approaches). You may alternatively provide a user-defined threshold value (lindley.thr option) to identify windows of interest\n")
    compute.localscore.pval=FALSE
  }else{
    compute.localscore.pval=TRUE
  }

  if(!is.null(lindley.thr)){
    if(lindley.thr<0){cat("WARNING: lindley.thr must be >0 (disregarded)\n")}
    if(compute.localscore.pval){
      cat("WARNING: the computation of p-value for the lindley-score will not be performed (the provided threshold will be used instead to identify windows)\n")
      compute.localscore.pval=FALSE}
  }

  ##Internal functions
  # computation of the autocorrelation
  autocor=function(x){abs(cor(x[-1],x[-length(x)]))}
  # computation of the lindley process from scores
  lindley=function(scores){
    L=length(scores) ; sl=rep(0,L+1)
    for (i in 1:L){sl[i+1]= max(0, (sl[i]+scores[i]))}
    return(sl[-1])
  }
  #computation of the significance threshold if the distribution of the p-values is uniform
  thresUnif=function(L, cor, xi, alpha = 0.05){
    a.poly.coef=rbind(c(-5.5,2.47,2.04,0.22),    #coeff associated with cor^3 for xi=1 to xi=4
                      c(6.76,-4.16,-5.76,-4.08), #coeff associated with cor^2 for xi=1 to xi=4
                      c(-5.66,-1.82,1.04,1.16), #coeff associated with cor   for xi=1 to xi=4
                      c(-2.51,-4.58,-6.95,-9.16))#intercept term              for xi=1 to xi=4
    b.poly.coef=rbind(c(-1.22,0.37,2.55,3.45),   #coeff associated with cor^2 for xi=1 to xi=4
                      c(3.17,2.14,-0.02,-0.98),  #coeff associated with cor   for xi=1 to xi=4
                      c(-1.99,-2.35,-2.31,-2.33)) #intercept term              for xi=1 to xi=4

    if(sum(xi != 1:4)==4){
      print('xi must be equal to 1, 2, 3 or 4')
      thres=NULL
    }else{#warning input cor might be a vector
      tmp.cor2=cor**2
      a=log(L)+a.poly.coef[1,xi]*(cor**3) + a.poly.coef[2,xi]*(cor**2) +
        a.poly.coef[3,xi]*cor + a.poly.coef[4,xi]
      b=b.poly.coef[1,xi]*(cor**2) + b.poly.coef[2,xi]*cor + b.poly.coef[3,xi]
      thres=( log(-log(1-alpha)) - a ) / b
    }
    return(thres)
  }
  ##

  if(is.bf){
    score=snp.bf[snp.sel]*bf.scale/10
    score[score<0]=-log10(runif(sum(score<0)))
  }else{
    score=snp.pvalue[snp.sel]
  }

  mydata=data.table(chr=snp.position[snp.sel,1],pos=snp.position[snp.sel,2],pval=10**(-score),score=score-xi)
  rm(score)
  setkey(mydata, chr)
  if(plot.pvalhist){
    hist(mydata$pval,breaks=100,freq=F,ylab="score",main="P-value distribution")
    abline(h=1,lty=3,col="red")
  }
  ##ordering of SNPs (in case not ordered) => pose pb (hard check plus haut: stop si un chr est pas ordonne)
  #  mydata <- mydata[order(chr, pos)]

  Nchr=length(mydata[,unique(chr)])
  #This is useful for doing genomewide plots (and compute thresholds)
  chrInf=mydata[,.(L=.N,cor=autocor(pval)),chr]
  setkey(chrInf,chr)
  # mydata[,score:=score-xi]
  # The score mean must be negative
  cat("Mean Score value:",mean(mydata$score),"\n")
  if(mean(mydata$score)>0){stop("The mean score must be negative (try increasing xi)\n")}
  mydata[,lindley:=lindley(score),chr]
  #significance threshold
  if(compute.localscore.pval){
    chrInf[,th01:=thresUnif(L, cor, xi=xi,alpha = pval.local.score.thres),]
  }else{
    if(!is.null(lindley.thr)){
      chrInf[,th01:=lindley.thr,]
      compute.localscore.pval=TRUE #pour pouvoir identifier les fenetres
    }
  }

  mydata=as.data.frame(mydata)

  if(manplot){
    tmp_chr=unique(mydata$chr)
    col_chr=1:length(tmp_chr) ; names(col_chr) = tmp_chr
    pos_chr = rep(0, length(tmp_chr))
    tmp_nmrk = table(mydata$chr)[tmp_chr] ;  pos_mrk = cumsum(tmp_nmrk)
    pos_chr[1] <- floor(pos_mrk[1]/2)
    if (length(tmp_chr) > 1) {
      for (i in 2:length(tmp_chr)) {pos_chr[i] = pos_mrk[i - 1] + floor((tmp_nmrk[i]/2))}
    }
    if(is.null(lindley.capping)){data=mydata$lindley}else{data=pmin(mydata$lindley,lindley.capping)}
    plot(data, pch = 16,las = 1,ylim=range(data,na.rm = TRUE),
         col = col_chr[mydata$chr], xaxt = "n",xlab = "",
         ylab = "Local-Score",...)
    axis(1, at = pos_chr, labels = tmp_chr, las = 1)
    chrInf=as.data.frame(chrInf)
    for(i in 1:nrow(chrInf)){lines(x=range(which(mydata$chr==chrInf$chr[i])),y=rep(chrInf$th01[i],2),col=i,lty=2)}
  }

  ###
  #detail region + annotation avec ma fonction (car autre ne sort pas la position du pic)
  #mydata=mydata[chrInf]
  #sigWin01=mydata[chrInf,sig_sl(lindley, pos, unique(th01)),chr]
  #ind=which(sigWin[,peak]>0)
  ls.win=matrix(NA,0,9)
  colnames(ls.win)=c("chr","beg","end","size","nsnps",
                     "Lindley score peak pos","Lindley score peak value",
                     "peak pos","-log10(p-val) peak value")
  if(is.bf){colnames(ls.win)[8:9]=c("BF (dB) peak pos","BF (dB) peak value")}
  if(compute.localscore.pval){
    for(i in 1:nrow(chrInf)){
      dum.sel=(mydata$chr==chrInf$chr[i])
      dum=mydata$lindley[dum.sel]
      dum.score.orig=mydata$score[dum.sel] + xi
      if(is.bf){dum.score.orig=snp.bf[snp.sel][dum.sel]}
      dum.pos=mydata[dum.sel,1:2]
      thr=chrInf$th01[i]
      while(max(dum)>thr){
        idx.max=which.max(dum)
        idx.init=rev(which(dum[1:idx.max]==0))[1]
        if(is.na(idx.init)){idx.init=1} #ca commence par >0
        idx.fin=idx.max+which(dum[-(1:idx.max)]==0)[1]
        if(is.na(idx.fin)){idx.fin=length(dum)}
        idx.max.score.orig=which.max(dum.score.orig[idx.init:idx.fin])
        dum.pos.beg=dum.pos$pos[idx.init]
        dum.pos.end=dum.pos$pos[idx.fin]
        dum.size=dum.pos$pos[idx.fin]-dum.pos$pos[idx.init]
        ls.win=rbind(ls.win,c(dum.pos$chr[idx.init],dum.pos.beg,dum.pos.end,dum.size,idx.fin-idx.init+1,
                              dum.pos$pos[idx.max],dum[idx.max],
                              dum.pos$pos[idx.init:idx.fin][idx.max.score.orig],
                              dum.score.orig[idx.init:idx.fin][idx.max.score.orig]))
        dum=dum[-(idx.init:idx.fin)]
        dum.score.orig=dum.score.orig[-(idx.init:idx.fin)]
        dum.pos=dum.pos[-(idx.init:idx.fin),]
      }
    }
    ls.win=as.data.frame(ls.win)
    for(i in 2:7){ls.win[,i]=as.numeric(ls.win[,i])}
  }
  return(list(res.local.scores=as.data.frame(mydata),significant.windows=as.data.frame(ls.win)))
}

