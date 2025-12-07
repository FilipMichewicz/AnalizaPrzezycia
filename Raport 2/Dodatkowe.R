# Wiktor Niedźwiedzki
# Filip Michewicz

library(coin)

times.II <- c(28, 89, 175, 195, 309, 377, 393, 421, 447, 462, 709, 744, 770, 1106, 1206)
deltas.II <- c(1,1,1,1,1,0,0,0,0,1,0,0,0,0,0)
df.II <- data.frame(times = times.II, deltas = deltas.II)

times.III <- c(34, 88, 137, 199, 280, 291, 299, 300, 309, 351, 358, 369, 369, 370, 375, 382, 392, 429, 451, 1119)
deltas.III <- c(1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1,1,0,1,0)
df.III <- data.frame(times = times.III, deltas = deltas.III)

group <- factor(c(rep("II", nrow(df.II)), rep("III", nrow(df.III))), levels = c("II", "III"))

df.all <- data.frame(times = c(df.II$times, df.III$times), deltas = c(df.II$deltas, df.III$deltas), group = group)

p_values_coin <- numeric(4)
names(p_values_coin) <- c("logrank", "Gehan-Breslow", "Tarone-Ware", "Peto-Peto")

for (name in names(p_values_coin)){
  p_values_coin[name] <- pvalue(logrank_test(Surv(times, deltas) ~ group, data = df.all, type = name))
}

### Testy od zera bo czemu by nie

event_times <- sort(unique(df.all$times[df.all$deltas == 1]))
groups <- unique(df.all$group)
k <- length(groups)

d_mat <- matrix(0, nrow = length(event_times), ncol = k)
r_mat <- matrix(0, nrow = length(event_times), ncol = k)

for(i in 1:length(event_times)){
  t <- event_times[i]
  for(g in 1:k){
    r_mat[i, g] <- sum(df.all$times >= t & df.all$group == groups[g])
    d_mat[i, g] <- sum(df.all$times == t & df.all$group == groups[g] & df.all$deltas == 1)
  }
}

r_sum <- rowSums(r_mat)
d_sum <- rowSums(d_mat)

S_minus <- numeric(length(event_times))
S_minus[1] <- 1
if(length(event_times) > 1){
  for(i in 2:length(event_times)){
    S_minus[i] <- S_minus[i-1] * (1 - d_sum[i-1]/r_sum[i-1])
  }
}

weights.df <- data.frame(
  logrank = rep(1, length(event_times)),
  Gehan_Breslow = r_sum,
  Tarone_Ware = sqrt(r_sum),
  Peto_Peto = S_minus
)

p_values_students <- numeric(4)
names(p_values_students) <- c("logrank", "Gehan_Breslow", "Tarone_Ware", "Peto_Peto")

for (name in names(weights.df)){
  W <- weights.df[[name]]
  Z <- numeric(k)
  for (j in 1:k){
    Z[j] <- sum(W * (d_mat[, j] - r_mat[, j] * d_sum / r_sum))
  }
  
  Z_reduced <- Z[1:(k-1)]
  sigma <- matrix(0, nrow = k-1, ncol = k-1)
  
  for(i in 1:(k-1)){
    for(j in 1:(k-1)){
      if(i == j){
        sigma[i,j] <- sum(W^2 * d_sum * r_mat[, j] / r_sum * (1 - r_mat[, j] / r_sum) * (r_sum - d_sum) / (r_sum - 1))
      } else {
        sigma[i,j] <- -sum(W^2 * d_sum * r_mat[, i] * r_mat[, j] / r_sum^2 * (r_sum - d_sum) / (r_sum - 1))
      }
    }
  }
  p_values_students[name] <- 1 - pchisq(t(Z_reduced) %*% solve(sigma) %*% Z_reduced, df = k-1)
}

# Porównanie

p_values_coin * 100
p_values_students * 100
