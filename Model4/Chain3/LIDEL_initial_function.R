LIDEL_initials<-function(mass, fs, fdoc, flig){
  
  #random draws for initial mass and fractions to C pools
  Initial_mass=rnorm(1, mean=mass[1], sd=mass[2])
  Initial_fs=rbeta(1, fs[1], fs[2])
  Initial_fdoc=rbeta(1, fdoc[1], fdoc[2])
  Initial_flig=rbeta(1, flig[1], flig[2])
  
  #calculate initial C pools
  LC1=Initial_fs*Initial_mass-(Initial_fs*Initial_fdoc*Initial_mass)
  LC2=(1-(Initial_fs+Initial_flig))*Initial_mass
  LC3=Initial_flig*Initial_mass
  LC6=Initial_fs*Initial_fdoc*Initial_mass
  
  #print(LC1+LC2+LC3+LC6)
  all=c(Initial_mass, LC1, LC2, LC3, LC6)
  names(all)=c("Mass", "C1", "C2", "C3", "C6")
  return(all)
}