# high spatial resolution differential analysis of RNA methylome with 
# fisher exact test and Hidden Markov Model
rhtestHMM <- function(untreated_ip, untreated_input, treated_ip, treated_input, 
                      untreated_ip_total, untreated_input_total, treated_ip_total, treated_input_total,
                      Pi=matrix(c(0.9,0.1,0.1,0.9),byrow=TRUE, nrow=2),
                      delta=c(0.5,0.5),
                      pm=list(prob=c(0.1, 0.9)),
                      threshold_fdr=0.05,
                      gene_label=rep(1,length(untreated_ip)),
                      minimal_count_fdr =10,mode="DIRECT"){
  # parameters
  # untreated_ip_total: an integer, sequencing depth of untreated ip sample   
  # untreated_input_total: ... untreated input sample                         
  # treated_ip_total: ... treated ip sample                                   
  # treated_input_total: ... treated input sample                             
  
  # untreated_ip:  a vector of integers, vector length equal to the number of binding sites, each element contains the number of reads within a binding site for the untreated ip sample
  # untreated_input:  ... of the untreated input sample
  # treated_ip:  ... of the treated ip sample
  # treated_input:  ... of the treated input sample
  
  # Pi:initial state transition probability matrix of HMM. default: matrix(c(0.9,0.1,0.1,0.9),byrow=TRUE, nrow=2)
  # delta:initial state probability of HMM. default: c(0.5,0.5)
  # pm:initial emission probability of HMM. default: c(0.1,0.9), i.e., the probability of observing a differential methylation state on an undifferential site is 0.1; and the probability for observing a differential methylation state on an differential methylation site is 0.9.
  # threshold_fdr:the threshold used to convert bltest states to observations of HMM model. default: 0.05, i.e., if a site is differential methylated with significance level 0.05, than it is considered a differential methylation loci in HMM model.
  # gene_label:if the counts are from multiple genes, then HMM is built on each gene. default: rep(1,length(untreated_ip)), i.e., the counts are all from the same gene.
  
  # minimal_count_fdr:an integer threshold, only the loci with reads more than this number are subjected for fdr calculation. default: 10
  # mode:a string,which specifies the mode of using HMM default:default:"DIRECT",which means to apply Estep of EM algorithm after estimating transition matrix and initial probability for all genes.
  # Alternative:"EM", which means to smooth the significance level(p-value) directly; "BEL",which means to use the Bernoulli process as the observation process of HMM.
  # required library
  require(exomePeak)
  require(HiddenMarkov)
  # rhtest analysis
  result<-rhtest(untreated_ip, untreated_input, treated_ip, treated_input, untreated_ip_total,
                 untreated_input_total, treated_ip_total, treated_input_total,
                 minimal_count_fdr = minimal_count_fdr)
  labelrhtest<- (result$log.fdr<log(threshold_fdr))
  if(mode=="DIRECT"){
    deltaestimate<-c(1-mean(labelrhtest),mean(labelrhtest))
    length<-length(treated_ip)
    state0<-c()
    p01<-c()
    p00<-c()
    state1<-c()
    p11<-c()
    p10<-c()
    # Pi
    state0 <- (labelrhtest[-length]==FALSE)
    p01<-mean(labelrhtest[-1][state0]) # 0 -> 1
    p00<-1- mean(labelrhtest[-1][state0]) # 0 -> 0
    state1 <- (labelrhtest[-length]==TRUE)
    p11<-mean(labelrhtest[-1][state1]) # 1 -> 1
    p10<-1 - mean(labelrhtest[-1][state1]) # 1 -> 0
    piestimate<-matrix(c(p00,p01,p10,p11), nrow = 2, ncol = 2, byrow = 2,dimnames = NULL)
    # HMM section
    pp<-exp(result$log.p)
    prob<-cbind(pp,1-pp)
    pen <- colSums(prob^2)
    prob[,1] <- prob[,1]/pen[1]
    prob[,2] <- prob[,2]/pen[2]
    unique_gl = unique(gene_label)
    pos <- cbind(untreated_ip*0 + 1,untreated_ip*0)
    #Estepresult<-list()
    for (i in 1:length(unique_gl)) {
      id = which(gene_label == unique_gl[i])
      m<-length(deltaestimate)
      len<-length(id)
      y<-.Estep_prob(piestimate, deltaestimate, m,len,prob[id,],fortran = TRUE, fwd.only = FALSE)
      pos[id,]<-y$u
    }
    # save HMM result
    m = untreated_ip + untreated_input + treated_ip + treated_input
    IDreads= which(m >minimal_count_fdr)
    #IDnot=which(m<= minimal_count_fdr)
    log.fdr<-c(rep(0,length(treated_ip)))
    log.p<-c(rep(0,length(treated_ip)))
    log.fdr[IDreads]=log(p.adjust(pos[IDreads,1], method = "fdr"))
    # remove infinity
    log.fdr=pmax(log.fdr,-1000)
    log.p=pmax(log(pos[,1]),-1000)
    #log.p[IDnot]=0
    # find peaks
    tmp <- ctest(IP=untreated_ip+treated_ip,INPUT=untreated_input+treated_input,
                 TOTAL_IP=untreated_ip_total+treated_ip_total,
                 TOTAL_INPUT=untreated_input_total+treated_input_total)
    tmp2 <- (tmp$log.fdr > log(0.05))
    # modify result
    log.fdr[tmp2] <- 0
    log.p[tmp2] <- 0
    result$log.fdr[tmp2]<-0
    result$log.p[tmp2]<-0
    fdr_fisher_sort <-sort(result$log.fdr)
    p_fisher_sort <-sort(result$log.p)
    log.fdr<- fdr_fisher_sort[rank(log.fdr)]
    log.p<- p_fisher_sort[rank(log.p)]
  }
  if(mode=="BEL"){
    # HMM section
    unique_gl = unique(gene_label)
    # initialization the posterior probability
    pos <- untreated_ip*0 + 1;
    for (i in 1:length(unique_gl)){
      id = which(gene_label == unique_gl[i])
      labelrhtest_i <- labelrhtest[id]
      pos[id] <- .HMM_single(labelrhtest_i,Pi,delta,pm)
    }
    # save HMM result
    log.fdr=log(p.adjust(pos, method = "fdr"))
    # adjust only sig reads count testing result
    m = untreated_ip + untreated_input + treated_ip + treated_input
    ID= which(m > minimal_count_fdr)
    log.fdr_sig=log(p.adjust(pos[ID], method = "fdr"))
    log.fdr[ID] = log.fdr_sig
    # remove infinity
    log.fdr=pmax(log.fdr,-1000)
    log.p=pmax(log(pos),-1000)
    
    }else if(mode=="EM"){
    # HMM section
    pp<-exp(result$log.p)
    prob<-cbind(pp,1-pp)
    # add the penalty to the p-values
    pen <- colSums(prob^2)
    prob[,1] <- prob[,1]/pen[1]
    prob[,2] <- prob[,2]/pen[2]
    unique_gl = unique(gene_label)
    # using the EM method(Baum-Welch Algorithm)
    rhtestHMM<-.Mstep_prob(Pi,delta,m,len, prob)
    # save HMM result
    log.fdr=log(p.adjust(rhtestHMM$gama, method = "fdr"))
    # adjust only sig reads count testing result
    m = untreated_ip + untreated_input + treated_ip + treated_input
    ID= which(m > minimal_count_fdr)
    log.fdr_sig=log(p.adjust(rhtestHMM$gama[ID], method = "fdr"))
    log.fdr[ID] = log.fdr_sig  
    # remove infinity
    log.fdr=pmax(log.fdr,-1000)
    log.p=pmax(log(rhtestHMM$gama),-1000)
  }
  # save result
  DIFF=list(log.fdr=log.fdr,log.p=log.p,log.fc=result$log.fc, rhtest_result = result)
  return(DIFF)
}
##########################################################################
# subfunction1
.HMM_single <- function(labelrhtest,Pi,delta,pm) {
  aa<-length(labelrhtest)
  pn <- list(size=rep(1,aa))
  x <- dthmm(labelrhtest, Pi, delta, "binom", pm, pn,discrete=TRUE)
  log <- capture.output({
    y <- BaumWelch(x);
  })
  pos <- y$u[,1]
  return(pos)
}
###############################################################################
#subfuction2
#define the Estep.prob function
.Estep_prob<-function(Pi, delta,m,len, pr, fortran = TRUE, fwd.only = FALSE){
  y <- forwardback.dthmm(Pi, delta, pr, fortran = TRUE, fwd.only = FALSE)
  logbeta <- y$logbeta
  logalpha <- y$logalpha
  LL <- y$LL
  u <- exp(logalpha + logbeta - LL)
  v <- array(NA, dim = c(len - 1, m, m))
  for (k in 1:m) {
    logPi <- matrix(log(Pi[, k]), byrow = TRUE, nrow = len -1, ncol = m)
    logPbeta <- matrix(log(pr[-1,k]) + logbeta[-1, k], byrow = FALSE, 
                       nrow = len - 1, ncol = m)
    v[, , k] <- logPi + logalpha[-len, ] + logPbeta - LL
  }
  v <- exp(v)
  return(list(u = u, v = v, LL = LL))
}
###############################################################################
#subfuction3(define the Mstep function)
.Mstep_prob<-function(pistart, deltastart,m,len, prob){
  pinew <- pistart
  deltanew <- deltastart
  iter<-1
  LL<-c()
  LL[1]=0
  # perform the E-step and M-step of Baum-Welch Algorithm
  repeat 
  { 
    m <- nrow(pinew)
    len<-length(untreated_ip)
    tsik<-c()
    tsij<-matrix(data = NA, nrow = 2, ncol = 2)
    y<-.Estep_prob(pinew, deltanew, m,len,prob, fortran = TRUE, fwd.only = FALSE)
    gama<- y$u
    ksi<- y$v
    LL[iter+1]<-y$LL
    # HMM M-step ,tisk & tsij are the expectation of times of transition.
    for (i in 1:2){
      tsik[i]=sum(gama[1:len-1,i])
      for (j in 1:2){
        tsij[i,j]=sum(ksi[1:len-1,i,j])}
    }
    # update delta and pi
    for (i in 1:2){
      deltanew[i]=gama[1,i]
      for (j in 1:2){
        pinew[i,j]=tsij[i,j]/tsik[i]}
    }
    print(iter)
    print(LL[iter])
    if(abs(LL[iter+1]-LL[iter])<(0.1^4))
    {
      break;
    }
    iter=iter+1;
  }
  #|| iter>100
  #the HMM converges
  iter<-iter
  gama<- y$u
  ksi<- y$v
  LLfinal<-y$LL
  return(list(gama = gama[,1], iter= iter, LLfinal = LLfinal))
}