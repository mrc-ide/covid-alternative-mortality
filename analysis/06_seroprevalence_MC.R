library(tidyverse)

#time to seroconversion
rexp(1,1/13.1)
hist(rexp(1000,1/13.1))
lines()


# time to seroreversion
rweibull(1,shape=3.7,scale=139)
hist(rweibull(1000,shape=3.7,scale=139))

# MC simulation

n <- 1000000

conv_time <- rexp(n, 1/13.1) #x (time between infection and conversion)
time_to_rev <- rweibull(n, shape = 3.7, scale = 139) #y (time between conversion and reversion)
rev_time <- conv_time + time_to_rev #x + y (time between infection and reversion)

hist(conv_time)
hist(time_to_rev,freq=FALSE)
hist(rev_time)

# construct individual timelines for each of the 100000 people so can determine each time when they would have antibodies or not

# time frame of our distribution (much longer than necessary)
upper_time <- 500

# test on one individual to see what happens before aggregating to n samples

# sample time to conversion
conv_time_ind <- round(rexp(1,1/13.1))
# record individual as negative until this time
pre_conv_time <- data.frame("time"=seq(0,conv_time_ind-1),"pos"=0)
# after being converting sample time to reversion
rev_time_ind <- round(rweibull(1, shape = 3.7, scale = 139))
# record individual as positive during this time
pos_time <- data.frame("time"=seq(conv_time_ind,rev_time_ind-1),"pos"=1)
# after reversion, individual is negative until the end of the sample time
neg_time_post_rev <- data.frame("time"=seq(rev_time_ind,upper_time),"pos"=0)
# combine data frames
status_ind <- rbind(pre_conv_time,pos_time,neg_time_post_rev)



ind_sim <- function(){
  conv_time_ind <- round(rexp(1,1/13.1))
  rev_time_ind <- round(rweibull(1, shape = 3.7, scale = 139))
  status_ind <- c(rep(0,times = length(seq(0,conv_time_ind-1))),
                  rep(1,times = length(seq(conv_time_ind,rev_time_ind-1))),
                  rep(0,times = length(seq(rev_time_ind,upper_time))))
  return(status_ind[1:501])
}

#ind_sim()

n_samp <- 100000

simulation <- data.frame(matrix(unlist(replicate(n=n_samp,ind_sim())),nrow=n_samp,byrow=TRUE))

plot(colSums(simulation)/n_samp,type="l",xlab="Time since infection",ylab="Probability of seropositivity",col="blue")
lines(colSums(simulation[1:100,])/100,type="l",col="red")
lines(colSums(simulation[1:1000,])/1000,type="l",col="green")
#lines(sero_det_igg,type="l",col="black") # this was the original function without the convolution (it matched)
legend(325, 0.3, legend=c("Simulation (100 samples)","Simulation (1000 samples)",
                          "Simulation (100000 samples)"#, "Original method"
                          ),
       col=c("red","green","blue"
             #, "black"
             ), lty=1:1, cex=0.8)


## plot for SI

pdfl <- 501

igg_conv <- 13.3; igm_conv <- 12.3
sero_conv_igg <- pexp(seq_len(pdfl), rate = 1/igg_conv)
sero_conv_igm <- pexp(seq_len(pdfl), rate = 1/igm_conv)

igg_scale <- 139; igm_scale <- 50

sero_rev_igg <- pweibull(seq_len(pdfl), 3.7, scale = igg_scale)
sero_rev_igm <- pweibull(seq_len(pdfl), 3.7, scale = igm_scale)

sero_det_igg <- sero_conv_igg - sero_rev_igg
sero_det_igg[sero_det_igg < 0] <- 0

sero_det_igm <- sero_conv_igm - sero_rev_igm
sero_det_igm[sero_det_igm < 0] <- 0



simulation_t <- t(simulation)

plot_sim <- data.frame("no"= seq(1,nrow(simulation_t),1),
                       "samp_100"=rowSums(simulation_t[,1:100])/100,
                       "samp_1000"=rowSums(simulation_t[,1:1000])/1000,
                       "samp_100000"=rowSums(simulation_t[])/ncol(simulation_t),
                       "original" = sero_det_igg)

ggplot(plot_sim %>% filter(no<=300))+
  geom_line(aes(x=no,y=original,col="Independent"),lwd=1.5)+
  geom_line(aes(x=no,y=samp_100,col="Dependent (100 samples)"),lwd=0.5)+
  geom_line(aes(x=no,y=samp_1000,col="Dependent (1000 samples)"),lwd=0.5)+
  geom_line(aes(x=no,y=samp_100000,col="Dependent (100000 samples)"),lwd=0.5)+
  theme_bw()+
  labs(x="Time since infection",y="Probability of seropositivity",col="")+
  scale_colour_manual(values = c("royalblue","#FF6666","#009900","black"))+
  theme(legend.position="bottom")+
  guides(col=guide_legend(nrow=2,byrow=TRUE))
ggsave("analysis/figures/seropositivity.png",height=4,width=6)


