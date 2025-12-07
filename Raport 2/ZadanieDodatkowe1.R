library(survminer)
library(survival)
library(coin)

times.II <- c(28, 89, 175, 195, 309, 377, 393, 421, 447, 462, 709, 744, 770, 
              1106, 1206)
deltas.II <- c(1,1,1,1,1,0,0,0,0,1,0,0,0,0,0)
groups.II <- rep(1, length(times.II))
df.II <- data.frame(times = times.II, deltas = deltas.II,
                    groups = groups.II)

times.III <- c(34, 88, 137, 199, 280, 291, 299, 300, 309, 351, 358, 369, 369, 
               370, 375, 382, 392, 429, 451, 1119)
deltas.III <- c(1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1,1,0,1,0)
groups.III <- rep(2, length(times.III))
df.III <- data.frame(times = times.III, deltas = deltas.III,
                     groups = groups.III)

list.df <- list(df.II[,-3], df.III[,-3])

df.cancer <- data.frame(rbind(df.II, df.III))
df.cancer$groups <- as.factor(df.cancer$groups)



test <- function(df_list, type="logrank"){
  
  # Łączenie prób w jedną
  df <- data.frame("times"=0, "deltas"=0)
  groups <- c()
  for(i in seq_along(df_list)){
    df <- rbind(df, df_list[[i]])
    groups <- c(groups, rep(i, nrow(df_list[[i]])))
  }

  df <- df[-1,]
  # Wyróżniamy unikalne czasy zdarzeń, dodajemy wartości r, d
  df.complete <- df[!duplicated(df$times) & df$deltas==1,]
  df.complete$r <- 0
  df.complete$d <- 0
  D <- nrow(df.complete)
  for(i in 1:D){
    df.complete$r[i] <- sum(df$times >= df.complete$times[i])
    df.complete$d[i] <- sum(df$times == df.complete$times[i] & df$deltas==1)
  }
  df$groups <- groups
  k <- length(unique(groups))
  for(g in 1:k){
    r_group <- numeric(D)
    d_group <- numeric(D)
    for(i in 1:D){
      r_group[i] <- sum(df$times >= df.complete$times[i] & df$groups == g)
      d_group[i] <- sum(df$times == df.complete$times[i] & 
                        df$deltas == 1 & df$groups == g)
    }
    col_name_r <- paste0("r", g)
    df.complete[[col_name_r]] <- r_group
    col_name_d <- paste0("d", g)
    df.complete[[col_name_d]] <- d_group
  }
  df.complete <- df.complete[order(df.complete$times),]

  # Liczenie wag
  if(type=="logrank") df.complete$weights <- 1
  else if(type=="Gehan-Breslo") df.complete$weights <- df.complete$r
  else if(type=="Tarone-Ware") df.complete$weights <- sqrt(df.complete$r)
  else if(type=="Peto-Peto"){
    weights <- numeric(D)
    product <- 1
    for(i in 1:(D-1)){
      weights[i] <- product
      product <- product * (1 - df.complete$d[i]/(df.complete$r[i]))
    }
    weights[D] <- product
    df.complete$weights <- weights
  }
  else if(type=="Prentice-Marek"){
    weights <- numeric(D)
    product <- 1
    for(i in 1:D){
      product <- product * (1 - df.complete$d[i]/(df.complete$r[i]+1))
      weights[i] <- product
    }
    df.complete$weights <- weights
  }
  else stop("Invalid option, choose one of the following: \"logrank\",
            \"Gehan-Breslo\", \"Tarone-Ware\", \"Peto-Peto\" or \"Prentice-Marek\".")
  
  # Wyznaczanie realizacji statystyk Z oraz estymatorów wariancji i kowariancji
  Sigma <- matrix(0, k-1, k-1)
  if(k-1 > 1){
    for(j in 1:(k-2)){
      r_j_column <- paste0("r", j)
      d_j_column <- paste0("d", j)
      for(g in (j+1):(k-1)){
        r_g_column <- paste0("r", g)
        d_g_column <- paste0("d", g)
        sum <- 0
        for(i in 1:D){
          r_i <- df.complete$r[i]
          difference <- r_i-1
          if(difference == 0) next 
          weight_i <- df.complete$weights[i]
          d_i <- df.complete$d[i]
          d_ij <- df.complete[[d_j_column]][i]
          r_ij <- df.complete[[r_j_column]][i]
          d_ig <- df.complete[[d_g_column]][i]
          r_ig <- df.complete[[r_g_column]][i]
          sum <- sum + weight_i^2*d_i*r_ij*r_ig/r_i^2*((r_i-d_i)/difference)
        }
        Sigma[j, g] <- -sum
      }
    }
  }
  Sigma <- Sigma + t(Sigma)
  Z <- numeric(k-1)
  variances <- numeric(k-1)
  for(j in 1:(k-1)){
    r_j_column <- paste0("r", j)
    d_j_column <- paste0("d", j)
    sum_Z <- 0
    sum_var <- 0
    for(i in 1:D){
      weight_i <- df.complete$weights[i]
      d_i <- df.complete$d[i]
      r_i <- df.complete$r[i]
      d_ij <- df.complete[[d_j_column]][i]
      r_ij <- df.complete[[r_j_column]][i]
      sum_Z <- sum_Z + weight_i * (d_ij - r_ij*d_i/r_i)
      difference <- r_i-1
      if(difference == 0) next
      sum_var <- sum_var +
                 weight_i^2*d_i*r_ij/r_i*(1-r_ij/r_i)*((r_i-d_i)/difference)
    }
    Z[j] <- sum_Z
    variances[j] <- sum_var
  }
  diag(Sigma) <- variances
  Z <- matrix(Z, nrow=1)
  
  # Obliczanie p-value
  p.value <- 1 - pchisq(as.numeric(Z %*% solve(Sigma) %*% t(Z)), df=k-1)
  return(p.value)
}


lg <- logrank_test(Surv(times, deltas) ~ groups,
                   data=df.cancer, type="logrank")
GB <- logrank_test(Surv(times, deltas) ~ groups,
                   data=df.cancer, type="Gehan-Breslo")
TW <- logrank_test(Surv(times, deltas) ~ groups,
                   data=df.cancer, type="Tarone-Ware")
PP <- logrank_test(Surv(times, deltas) ~ groups,
                   data=df.cancer, type="Peto-Peto")
PM <- logrank_test(Surv(times, deltas) ~ groups,
                   data=df.cancer, type="Prentice-Marek")

test_lg <- test(list.df, type="logrank")
test_GB <- test(list.df, type="Gehan-Breslo")
test_TW <- test(list.df, type="Tarone-Ware")
test_PP <- test(list.df, type="Peto-Peto")
test_PM <- test(list.df, type="Prentice-Marek")

# logrank
pvalue(lg)
test_lg
test(rev(list.df), type="logrank")
abs(pvalue(lg)-test_lg)

# Gehan-Breslow
pvalue(GB)
test_GB
test(rev(list.df), type="Gehan-Breslo")
abs(pvalue(GB)-test_GB)

# Tarone-Ware
pvalue(TW)
test_TW
test(rev(list.df), type="Tarone-Ware")
abs(pvalue(TW)-test_TW)

# Peto-Peto
pvalue(PP)
test_PP
test(rev(list.df), type="Peto-Peto")
abs(pvalue(PP)-test_PP)

# Prentice-Marek
pvalue(PM)
test_PM
test(rev(list.df), type="Prentice-Marek")
abs(pvalue(PM)-test_PM)



max(abs(pvalue(lg)-test_lg),
    abs(pvalue(GB)-test_GB),
    abs(pvalue(TW)-test_TW),
    abs(pvalue(PP)-test_PP),
    abs(pvalue(PM)-test_PM))