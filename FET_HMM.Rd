\name{rhtestHMM}
\alias{rhtestHMM}
\title{
rhtestHMM
}
\description{
This function combines rhtest and HMM, which aims at high spatial resolution differential analysis of RNA methylome.
}
\usage{
rhtestHMM(untreated_ip, untreated_input, 
       treated_ip, treated_input, 
       untreated_ip_total, untreated_input_total, 
       treated_ip_total, treated_input_total,
       Pi=matrix(c(0.9,0.1,0.1,0.9),byrow=TRUE, nrow=2),
       delta=c(0.5,0.5),pm=list(prob=c(0.1, 0.9)),
       threshold_fdr=0.05,gene_label=rep(1,length(untreated_ip)),
       minimal_count_fdr =10,mode="BEL")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{untreated_ip}{ a vector of integers of n, which is the number of binding sites tested. Each element represents the number of reads fall into a binding site for the IP sample under untreated condition
%%     ~~Describe \code{untreated_ip} here~~
}
  \item{untreated_input}{ a vector of integers of n, which is the number of binding sites tested. Each element represents the number of reads fall into a binding site for the Input control sample under untreated condition
%%     ~~Describe \code{untreated_input} here~~
}
  \item{treated_ip}{ a vector of integers of n, which is the number of binding sites tested. Each element represents the number of reads fall into a binding site for the IP sample under treated condition
%%     ~~Describe \code{treated_ip} here~~
}
  \item{treated_input}{ a vector of integers of n, which is the number of binding sites tested. Each element represents the number of reads fall into a binding site for the Input control sample under treated condition
%%     ~~Describe \code{treated_input} here~~
}
  \item{untreated_ip_total}{ an integer, total number of reads for the IP sample under untreated condition
%%     ~~Describe \code{untreated_ip_total} here~~
}
  \item{untreated_input_total}{ an integer, total number of reads for the Input control sample under untreated condition
%%     ~~Describe \code{untreated_input_total} here~~
}
  \item{treated_ip_total}{ an integer, total number of reads for the IP sample under treated condition
%%     ~~Describe \code{treated_ip_total} here~~
}
  \item{treated_input_total}{  an integer, total number of reads for the Input control sample under treated condition
%%     ~~Describe \code{treated_input_total} here~~
}
\item{Pi}{ initial state transition probability matrix of HMM. default: matrix(c(0.9,0.1,0.1,0.9),byrow=TRUE, nrow=2)
%%     ~~Describe \code{treated_input_total} here~~
}
\item{delta}{ initial state probability of HMM. default: c(0.5,0.5).
%%     ~~Describe \code{treated_input_total} here~~
}
\item{pm}{initial emission probability of HMM. default: c(0.1,0.9), i.e., the probability of observing a differential methylation state on an undifferential site is 0.1; and the probability for observing a differential methylation state on an differential methylation site is 0.9.
%%     ~~Describe \code{treated_input_total} here~~
}
\item{threshold_fdr}{ the threshold used to convert bltest states to observations of HMM model. default: 0.05, i.e., if a site is differential methylated with significance level 
0.05, than it is considered a differential methylation loci in HMM model.
%%     ~~Describe \code{treated_input_total} here~~
}
\item{gene_label}{ if the counts are from multiple genes, then HMM is built on each gene. default: rep(1,length(untreated_ip)), i.e., the counts are all from the same gene.
%%     ~~Describe \code{treated_input_total} here~~
}
 \item{minimal_count_fdr}{ an integer threshold, only the loci with reads more than this number are subjected for fdr calculation. default: 10
%%     ~~Describe \code{treated_input_total} here~~
}
 \item{mode}{ a string,which specifies the mode of using HMM. default:"DIRECT",which means to apply Estep of EM algorithm after estimating transition matrix and initial probability for all genes.Alternative:"BEL",which means to use the Bernoulli process as the observation process of HMM;"EM",which means to smooth the significance level(p-value) directly.
%%     ~~Describe \code{treated_input_total} here~~
}
}
\value{
The function returns a list of length 4, which contains the log(p-value), log(fdr) ,log(fold change)and the initial result of "rhtest",respectively, from the analysis.
}

\details{
This function combines rhtest and HMM.Based on the result of ¡®RHtest¡¯ analysis, we applied HMM to conduct the further high resolution analyzing step.The number of reads at the same locus for the aligned reads are counted by other packages, such as Rsamtools or HTseq-count.
}

\references{
Reference coming soon!
}
\author{ Yuchen Zhang, M.S., <yuchenzhang@mail.nwpu.edu.cn>
%%  ~~who you are~~
}

\examples{
# these 4 vectors are generated to simulate MeRIP-Seq data with 4 adjacent sections, 
# including 2 undifferential sites (US), and 2 differential methylation sites (DMS).
# the order of them is: US-DMS-US-DMS
untreated_ip = c(rpois(30,20),rpois(30,60),rpois(30,20), rpois(30,5)) 
untreated_input = rep(20,120)
treated_ip = rep(20,120)
treated_input = rep(20,120)
# sequencing depths
untreated_ip_total = 10^7
untreated_input_total = 10^7
treated_ip_total = 10^7
treated_input_total = 10^7
result <- rhtestHMM(untreated_ip, untreated_input, 
                    treated_ip, treated_input, 
                    untreated_ip_total, untreated_input_total, 
                    treated_ip_total, treated_input_total,mode="DIRECT")
# prediction
prediction <- (rbind(hmm=result$log.fdr,
                     rhtest=result$rhtest_result$log.fdr) < log(0.05))
# comparison
comparison<- rbind(setting = c(rep(FALSE,30),rep(TRUE,30),rep(FALSE,30),rep(TRUE,30)), 
             prediction);
hmm_perf <- 1-mean(abs(comparison[1,]-comparison[2,]))
rhtest_perf <- 1-mean(abs(comparison[1,]-comparison[3,]))
compare_result <- data.frame(rhtestHMM_performance = hmm_perf, rhtest_performance=rhtest_perf)
# show the performance of two methods
print(compare_result)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ bioinformatics statistical analysis}
