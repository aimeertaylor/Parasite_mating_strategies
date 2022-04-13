#' #################################################################################
#' Simulation of exclusive superinfection (prohibits brood mating in next
#' mosquito) versus exclusive cotransmission (prohibits non-brood mating in the
#' next mosquito) versus mixed superinfection and cotransmission as a function
#' of number of sporozoites innoculated and various other parameters.
#' 
#' Make shiny
#' Check monoclonal
#' Check MOI increases with bites
#' Remember we don't have data on people who were bitten but didn't develop an infection
#' Check if probability of monoclonal infection is biased or not
#' #################################################################################
rm(list = ls())

#' Function that simulates the number of genotypes in a single bite per person
bite_per_person <- function(n_p, n_i, p_mono, p_bld){
  
  #' Vector of MOIs
  n_geno <- rep(NA, n_p)
  
  #' Index of people bitten by monolclonally infected
  #' mosquitoes among the n_p people bitten
  i_mono <- as.logical(rbinom(n = n_p, size = 1, prob = p_mono))
  
  #' Number of sporozoites successivly resulting in a blood-stage infection for 1
  #' to m_mono mosquito bites
  n_mono <- rbinom(n = sum(i_mono == T), size = n_i, prob = p_bld) 
  
  # Record number of genotypes for people bitten by monocloanlly infected mosquitoes
  n_geno[i_mono] <- 1*(n_mono > 0)
  
  # If there are any multiply infected mosquitoes... 
  if(sum(i_mono == F) > 0){
    
    #' Number of sporozoites successly resulting in a blood-stage infection for 1
    #' to m_poly mosquito bites
    n_poly <- rbinom(n = sum(i_mono == F), size = n_i, prob = p_bld)
    
    #' Draw genotypes for n_poly sporozoites in 1 to m_poly mosquito bites with
    #' zero or more successful sporozoite
    g_poly <- sapply(n_poly, function(x){
      rmultinom(n = 1, size = x, prob = rep(1/(4*n_o), (4*n_o)))})
    
    #' Record number of genotypes for people bitten by monocloanlly infected mosquitoes
    n_geno[!i_mono] <- apply(g_poly, 2, function(x) sum(x > 0)) 
  }
  
  if(any(is.na(n_geno))) {stop("something wrong")}
  return(n_geno)    
}


# Among twice bitten individuals
sim <- function(n_bpp, n_p, n_is, p_mono, p_bld, n_o){
  
  simulation <- sapply(n_is, function(n_i){
    
    # number of distinct genotypes following the 1:n_bpp bites per n_p people 
    X <- as.matrix(sapply(1:n_bpp, FUN = function(x){bite_per_person(n_p, n_i, p_mono, p_bld)}))
    
    # logical for polyclonal / monoclonal
    mois <- rowSums(X)
    polyclonal <- mois > 1
    monoclonal <- mois == 1
    
    # Proportion of people with monoclonal infection 
    # and average moi 
    # among those with with a blood stage infection
    # (for comparison with observed data)
    prop_mono <- sum(monoclonal)/(sum(polyclonal)+sum(monoclonal))
    moi_av <- mean(mois[mois > 0])
    
    # Proportion MOI equal to or less than number of bites per person 
    moi_nbpp_equal <- mean(mois <= n_bpp)
    
    # fully co-transmitted
    co <- mean(rowSums(X[polyclonal,, drop = F] != 0) == 1)
    
    # fully superinfection
    sup <- mean(apply(X[polyclonal,, drop = F], 1, function(x) all(x <= 1)))
    
    c(Cotransmission = co, 
      Superinfection = sup, 
      MOI_nbpp = moi_nbpp_equal, # Multiplicity of infection <= nbpp
      MOI_max_pb = max(X), # Max MOI per bite 
      MOI_av = moi_av, # Multiplicity of infection average
      prop_mono = prop_mono) # Proportion of monoclonally infected among infected
  })
  
  return(simulation)
}


X <- sapply(1:5, simplify = "array", FUN = function(n_bpp){
  n_is = seq(1,1000,10) # Numbers of sporozoites innoculated
  n_p = 5000 # Number of people who receive a bite from a mosquito who innoculate n_i sporozoites (number of replicates)
  p_mono = 0.3 # Probability of monoclonally infected mosquito
  p_bld = 0.02 # Probability of sporozoite resulting in a blood-stage infection
  n_o = 2 # Number of ooycsts (little effect, presumably because 2 is enough)
  sim(n_bpp, n_p, n_is, p_mono, p_bld, n_o)
  })

# Check that with EIR 1) monoclonal decreases and 2) average MOI increases
moi_max <- ceiling(max(X["MOI_av",,]))
par(mar = c(4,4,2,5))
plot(NULL, ylim = c(0,1), xlim = range(n_is), 
     xlab = "Number of innocualted sporozoites", 
     ylab = "Proportion of monoclonal blood-stage infections")
mtext(side = 4, "MOI", las = 2, line = 3)
rect(xleft = 10, xright = 100, ytop = 1, ybottom = 0, col = "lightgray", border = NA)
axis(side = 4, labels = seq(1:moi_max), at = ((1:moi_max)-1)/(moi_max - 1), las = 1)
for(n_bpp in 1:5){
  print(max(X["MOI_av",,n_bpp]))
  lines(y = (X["MOI_av",,n_bpp]-1)/(moi_max - 1), x = n_is, col = n_bpp)
  lines(y = X["prop_mono",,n_bpp], x = n_is, col = n_bpp, lty = "dashed")
}
legend("top", bty = 'n', legend = 1:5, title = "Number bites per person", 
       fill = 1:5, cex = 0.5, inset = 0.03, border = NA)
legend("top", bty = "n", lty = c("dashed","solid"), 
       legend = c("Proportion monoclonal","
                  MOI"), 
       cex = 0.5, inset = 0.15)


plot(NULL, ylim = c(0,1), xlim = c(10,100), 
     xlab = "Number of innocualted sporozoites", 
     ylab = "Proportion of monoclonal blood-stage infections")
mtext(side = 4, "MOI", las = 2, line = 3)
rect(xleft = 10, xright = 100, ytop = 1, ybottom = 0, col = "lightgray", border = NA)
axis(side = 4, labels = seq(1:moi_max), at = ((1:moi_max)-1)/(moi_max - 1), las = 1)
for(n_bpp in 1:5){
  print(max(X["MOI_av",,n_bpp]))
  lines(y = (X["MOI_av",,n_bpp]-1)/(moi_max - 1), x = n_is, col = n_bpp)
  lines(y = X["prop_mono",,n_bpp], x = n_is, col = n_bpp, lty = "dashed")
}
legend("top", bty = 'n', legend = 1:5, title = "Number bites per person", 
       fill = 1:5, cex = 0.5, inset = 0.03, border = NA)
legend("top", bty = "n", lty = c("dashed","solid"), 
       legend = c("Proportion monoclonal","MOI"), 
       cex = 0.5, inset = 0.15)

par(mar = c(4,6,2,1))
plot(NULL, xlim = range(10,100), ylim = c(0,1.1), bty = "n", panel.first = grid(),
     ylab = sprintf("Among %s times innoculated individuals,
                    proportion of polyclonal infections whose parasites were
                    exclusively co-transmitted or exclusively superinfected or neither", n_bpp), 
     xlab = "Number of innoculated sporozoites", cex.lab = 0.75)
rect(xleft = 10, xright = 100, ytop = 1, ybottom = 0, col = "gray", border = NA)
lines(y = X["Cotransmission",,2], x = n_is, col = "red")
lines(y = X["Superinfection",,2], x = n_is, col = "blue")
lines(y = 1-colSums(X[c("Cotransmission", "Superinfection"),,2]), 
      x = n_is, col = "purple")
lines(y = X["prop_mono",], x = n_is, col = "black")
legend('right', 
       bty = 'n', legend = c("Exclusively co-transmitted", 
                             "Exclusively superinfected", 
                             "Mixture"), cex = 0.5, 
       col = c("red", "blue", "purple"), lty = 1)


#====================================================================

for(n_bpp in 1:5){
  
  p_mono <- 0 # Probability of monoclonally infected mosquito
  p_bld <- 0.01 # Probability of sporozoite resulting in a blood-stage infection
  n_o <- 1 # Number of ooycsts (little effect, presumably because 2 is enough)
  
  X <- sim(n_bpp, n_p, n_is, p_mono, p_bld, n_o)
  
  png(filename = sprintf("~/Desktop/sim_pbld%s_pmono%s_no%s_nbpp%s.png", p_bld, p_mono, n_o, n_bpp), 
      width = 1000, height = 1000, res = 200)
  par(mar = c(4,6,3,1))
  plot(NULL, xlim = range(n_is), ylim = c(0,1.1), bty = "n", panel.first = grid(),
       ylab = sprintf("Among %s times innoculated individuals,
                      proportion of polyclonal infections whose parasites were
                      exclusively co-transmitted or exclusively superinfected or neither", n_bpp), 
       xlab = "Number of innoculated sporozoites", cex.lab = 0.75)
  rect(xleft = 10, xright = 100, ytop = 1, ybottom = 0, col = "gray", border = NA)
  lines(y = X["Cotransmission",], x = n_is, col = "red")
  lines(y = X["Superinfection",], x = n_is, col = "blue")
  lines(y = 1-colSums(X[c("Cotransmission", "Superinfection"),]), 
        x = n_is, col = "purple")
  lines(y = X["prop_mono",], x = n_is, col = "black")
  legend('right', 
         bty = 'n', legend = c("Exclusively co-transmitted", 
                               "Exclusively superinfected", 
                               "Mixture"), cex = 0.5, 
         col = c("red", "blue", "purple"), lty = 1)
  
  mtext(side = 3, line = -1, cex = 0.5, text = sprintf("
                                                       Probability of sporozoite resulting in a blood-stage infection %s, 
                                                       Probability of being bitten by a monoclonally infected mosquito %s, 
                                                       Number of infective bites per person over infection period %s, 
                                                       Number of ooycsts %s", p_bld, p_mono, n_bpp, n_o))
  dev.off()
  
  plot(y = X["MOI_av",], x = n_is, col = "black")
}



for(n_o in 1:5){
  
  p_mono <- 0 # Probability of monoclonally infected mosquito
  p_bld <- 0.01 # Probability of sporozoite resulting in a blood-stage infection
  n_bpp <- 2 # Number of infective bites per person over untreated/unrecovered period
  
  X <- sim(n_bpp, n_p, n_is, p_mono, p_bld, n_o)
  
  png(filename = sprintf("~/Desktop/sim_pbld%s_pmono%s_no%s_nbpp%s.png", p_bld, p_mono, n_o, n_bpp), 
      width = 1000, height = 1000, res = 200)
  par(mar = c(4,6,3,1))
  plot(NULL, xlim = range(n_is), ylim = c(0,1.1), bty = "n", panel.first = grid(),
       ylab = sprintf("Among %s times innoculated individuals,
                        proportion of polyclonal infections whose parasites were
                        exclusively co-transmitted or exclusively superinfected or neither", n_bpp), 
       xlab = "Number of innoculated sporozoites", cex.lab = 0.75)
  rect(xleft = 10, xright = 100, ytop = 1, ybottom = 0, col = "gray", border = NA)
  lines(y = X["Cotransmission",], x = n_is, col = "red")
  lines(y = X["Superinfection",], x = n_is, col = "blue")
  lines(y = 1-colSums(X[c("Cotransmission", "Superinfection"),]), 
        x = n_is, col = "purple")
  legend('right', 
         bty = 'n', legend = c("Exclusively co-transmitted", 
                               "Exclusively superinfected", 
                               "Mixture"), cex = 0.5, 
         col = c("red", "blue", "purple"), lty = 1)
  
  mtext(side = 3, line = -1, cex = 0.5, text = sprintf("
                                                         Probability of sporozoite resulting in a blood-stage infection %s, 
                                                         Probability of being bitten by a monoclonally infected mosquito %s, 
                                                         Number of infective bites per person over infection period %s, 
                                                         Number of ooycsts %s", p_bld, p_mono, n_bpp, n_o))
  dev.off()
}

for(p_mono in seq(0, 1, 0.2)){
  
  p_bld <- 0.01 # Probability of sporozoite resulting in a blood-stage infection
  n_o <- 1 # Number of ooycsts (little effect, presumably because 2 is enough)
  n_bpp <- 2 # Number of infective bites per person over untreated/unrecovered period
  
  X <- sim(n_bpp, n_p, n_is, p_mono, p_bld, n_o)
  
  png(filename = sprintf("~/Desktop/sim_pbld%s_pmono%s_no%s_nbpp%s.png", p_bld, p_mono, n_o, n_bpp), 
      width = 1000, height = 1000, res = 200)
  par(mar = c(4,6,3,1))
  plot(NULL, xlim = range(n_is), ylim = c(0,1.1), bty = "n", panel.first = grid(),
       ylab = sprintf("Among %s times innoculated individuals,
                      proportion of polyclonal infections whose parasites were
                      exclusively co-transmitted or exclusively superinfected or neither", n_bpp), 
       xlab = "Number of innoculated sporozoites", cex.lab = 0.75)
  rect(xleft = 10, xright = 100, ytop = 1, ybottom = 0, col = "gray", border = NA)
  lines(y = X["Cotransmission",], x = n_is, col = "red")
  lines(y = X["Superinfection",], x = n_is, col = "blue")
  lines(y = 1-colSums(X[c("Cotransmission", "Superinfection"),]), 
        x = n_is, col = "purple")
  legend('right', 
         bty = 'n', legend = c("Exclusively co-transmitted", 
                               "Exclusively superinfected", 
                               "Mixture"), cex = 0.5, 
         col = c("red", "blue", "purple"), lty = 1)
  
  mtext(side = 3, line = -1, cex = 0.5, text = sprintf("
                                                       Probability of sporozoite resulting in a blood-stage infection %s, 
                                                       Probability of being bitten by a monoclonally infected mosquito %s, 
                                                       Number of infective bites per person over infection period %s, 
                                                       Number of ooycsts %s", p_bld, p_mono, n_bpp, n_o))
  dev.off()
}


for(p_bld in rev(c(0.001, 0.005, 0.01, 0.05, 0.1))){
  
  p_mono <- 0 # Probability of monoclonally infected mosquito
  n_o <- 1 # Number of ooycsts (little effect, presumably because 2 is enough)
  n_bpp <- 2 # Number of infective bites per person over untreated/unrecovered period
  
  X <- sim(n_bpp, n_p, n_is, p_mono, p_bld, n_o)
  
  png(filename = sprintf("~/Desktop/sim_pbld%s_pmono%s_no%s_nbpp%s.png", p_bld, p_mono, n_o, n_bpp), 
      width = 1000, height = 1000, res = 200)
  par(mar = c(4,6,3,1))
  plot(NULL, xlim = range(n_is), ylim = c(0,1.1), bty = "n", panel.first = grid(),
       ylab = sprintf("Among %s times innoculated individuals,
                      proportion of polyclonal infections whose parasites were
                      exclusively co-transmitted or exclusively superinfected or neither", n_bpp), 
       xlab = "Number of innoculated sporozoites", cex.lab = 0.75)
  rect(xleft = 10, xright = 100, ytop = 1, ybottom = 0, col = "gray", border = NA)
  lines(y = X["Cotransmission",], x = n_is, col = "red")
  lines(y = X["Superinfection",], x = n_is, col = "blue")
  lines(y = 1-colSums(X[c("Cotransmission", "Superinfection"),]), 
        x = n_is, col = "purple")
  legend('right', 
         bty = 'n', legend = c("Exclusively co-transmitted", 
                               "Exclusively superinfected", 
                               "Mixture"), cex = 0.5, 
         col = c("red", "blue", "purple"), lty = 1)
  
  mtext(side = 3, line = -1, cex = 0.5, text = sprintf("
Probability of sporozoite resulting in a blood-stage infection %s, 
Probability of being bitten by a monoclonally infected mosquito %s, 
Number of infective bites per person over infection period %s, 
Number of ooycsts %s", p_bld, p_mono, n_bpp, n_o))
  dev.off()
  
}



# # Checks 
# # MOI_nbpp is an upper bound on superinfection
# lines(y = sim["MOI_nbpp",], x = n_is, col = "blue", lty = "dashed")
# # MOIno_max is the maxium muliplicity of innoculation - check it increases with n_o
# lines(y = (sim["MOI_max_pb",]-min(sim["MOI_max_pb",]))/diff(range(sim["MOI_max_pb",])), x = n_is, col = "black" )




