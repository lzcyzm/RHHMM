# example of using bltestHMM
untreated_ip = c(rpois(30,20),rpois(30,60),rpois(30,20), rpois(30,5)) 
# this is to simulate MeRIP-Seq data with 4 adjacent sections, 
# including 2 undifferential sites (US), and 2 differential methylation sites (DMS).
# the order of them is: US-DMS-US-DMS
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