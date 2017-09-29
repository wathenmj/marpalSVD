# Rao_canonical_phi
a<-0.3;b<-0.4;c<-0.2;d<-0.1
n<-1

MB <- rdirichlet(1, c(1,1,1,1))
Sample<-rmultinom(1, 100, prob= c(MB[1,1], MB[1,2], MB[1,3], MB[1,4]))


# assymptotic variance of phi_hat

variance_phi_hat<-function(a,b,c,d){
  population_data<-matrix(c(a,b,c,d),nrow=2,ncol=2,byrow=T)

  population_gamma<-(a+b)*(a+c)*(c+d)*(b+d)
  phi<-((a*d) - (b*c)) * (-sqrt(population_gamma))

  v_A<- a * (1-a)
  v_B<- b * (1-b)
  v_C<- c * (1-c)
  v_D<- d * (1-d)
  cov_ab<- a * b
  cov_ac<- a * c
  cov_ad<- a * d
  cov_bc<- b * c
  cov_bd<- b * d
  cov_cd<- c * d

  df_da<-((1/n)*d*(-sqrt(population_gamma))) - (0.5*phi*((1+a+d)/((a+b)*(a+c))))
  df_db<-((-1/n)*c*(-sqrt(population_gamma))) + (0.5*phi*((1+b-c)/((b+d)*(b+a))))
  df_dc<-((-1/n)*b*(-sqrt(population_gamma))) + (0.5*phi*((1+c-b)/((c+a)*(c+d))))
  df_dd<-((1/n)*a*(-sqrt(population_gamma))) - (0.5*phi*((1+d-a)/((d+b)*(d+c))))

  V_phi_hat <- v_A * (df_da)^2 +  v_B * (df_db)^2 + v_C * (df_dc)^2 +
    v_D * (df_dd)^2 + 2 * cov_ab * (df_da) * (df_db) +
    2 * cov_ac * (df_da) * (df_dc) + 2 * cov_ad * (df_da) * (df_dd)+
    2 * cov_bc * (df_db) * (df_dc) + 2 * cov_bd * (df_db) * (df_dd)+
    2 * cov_cd * (df_dc) * (df_dd)


  return (list(Data=population_data,Variance=V_phi_hat))

}

