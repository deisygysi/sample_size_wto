## Characterize the convergence for the pvalue for
## the wto

require(wTO)
require(data.table)
require(magrittr)
set.seed(123)

### wTO function
wTO = function(A){
  A_TF = A
  C = as.matrix(A) %*% t(A)
  W = C + A_TF
  K = matrix(NA, nrow(A_TF), ncol(A_TF))
  KI = rowSums(abs(A), na.rm = T)
  for (ii in 1:nrow(A_TF)) {
    for (jj in 1:ncol(A_TF)) {
      K[ii, jj] = min(KI[ii], KI[jj])
    }
  }
  WTO = round(W/(K + 1 - abs(A_TF)), 3)
  return(WTO)
}

mat = rnorm(1000) %>%
matrix(ncol = 10, nrow = 10) %>%
as.data.frame()
row.names(mat)= LETTERS[1:10]

mat = wTO::Microarray_Expression1[10:100,]
Data = mat

A = wTO::CorrelationOverlap(Data = Data, Overlap = row.names(Data),
                            method = "p") %>% as.data.frame()

N = 10
O = list()
O[[1]] = wTO(A) %>% wTO.in.line()
for ( i in 1:N){
  A = wTO::CorrelationOverlap(Data =  Data[, sample(1:ncol(Data),
                                                    replace = TRUE)], Overlap = row.names(Data),
                              method = "p") %>% as.data.frame()
  O[[i+1]] = wTO(A) %>% wTO.in.line()
  names(O[[i+1]])[3]=paste0("C_", i)
}

C = O%>% plyr::join_all()
D = C[, -c(1:3)]
D[1,] %>% as.numeric() %>% hist
abline(v = C[1,3], col = "red")

x_10 = C[1, ]
y_10 = C[92,]
###
N = 100
O = list()
O[[1]] = wTO(A) %>% wTO.in.line()
for ( i in 1:N){
  A = wTO::CorrelationOverlap(Data =  Data[, sample(1:ncol(Data),
                                                    replace = TRUE)], Overlap = row.names(Data),
                              method = "p") %>% as.data.frame()
  O[[i+1]] = wTO(A) %>% wTO.in.line()
  names(O[[i+1]])[3]=paste0("C_", i)
}

C = O%>% plyr::join_all()
D = C[, -c(1:3)]
D[1,] %>% as.numeric() %>% hist
abline(v = C[1,3], col = "red")
x_100 = C[1,]
y_100 = C[92,]

N = 1000
O = list()
O[[1]] = wTO(A) %>% wTO.in.line()
for ( i in 1:N){
  A = wTO::CorrelationOverlap(Data =  Data[, sample(1:ncol(Data),
                                                    replace = TRUE)], Overlap = row.names(Data),
                              method = "p") %>% as.data.frame()
  O[[i+1]] = wTO(A) %>% wTO.in.line()
  names(O[[i+1]])[3]=paste0("C_", i)
}

### Select random links to check convergence
C = O%>% plyr::join_all()
D = C[, -c(1:3)]
D[1,] %>% as.numeric() %>% hist
abline(v = C[1,3], col = "red")
x_1000 = C[1,]
y_1000 = C[92,]

a = x_10[,-c(1:3)] %>% t() %>% as.data.frame()
b = x_100[,-c(1:3)] %>% t() %>% as.data.frame()
c = x_1000[,-c(1:3)] %>% t() %>% as.data.frame()

a$sample = 10
b$sample = 100
c$sample = 1000

plot_ = rbind(a,b,c)
plot_$sample %<>% as.factor()

x = C[1,3] %>% as.numeric()
require(ggplot2)
ggplot(plot_) +
  aes(x = V1, fill = sample, colour = sample) +
  geom_density(alpha = 0.1) +
  scale_fill_hue() +
  scale_color_hue() +
  geom_vline(xintercept = x)+
  theme_minimal()


example_2 =C[1, -c(1:3)] %>% as.numeric() %>% data.frame()
names(example_2) = "wTO"
example_2$run = 1:1000

delta = 0.2
min = (C[1,3] - delta) %>% as.numeric()
max = (C[1,3] + delta) %>% as.numeric()
example_2$p = ifelse( example_2$wTO > min & example_2$wTO < max, 0, 1)

example_2$prop = cumsum(example_2$p)
example_2$pval = example_2$prop/example_2$run


ggplot(example_2) +
  aes(x = run, y = pval) +
  geom_point(size = 1L, colour = "#bd3786") +
  geom_vline(xintercept = 100)+
  # geom_smooth(span = 1L, colour = "#bd3786") +
  theme_minimal()




#####
d = y_10[,-c(1:3)] %>% t() %>% as.data.frame()
e = y_100[,-c(1:3)] %>% t() %>% as.data.frame()
f = y_1000[,-c(1:3)] %>% t() %>% as.data.frame()

d$sample = 10
e$sample = 100
f$sample = 1000

plot_ = rbind(d,e,f)
plot_$sample %<>% as.factor()

x = C[92,3] %>% as.numeric()
require(ggplot2)
ggplot(plot_) +
  aes(x = V1, fill = sample, colour = sample) +
  geom_density(alpha = 0.1) +
  scale_fill_hue() +
  scale_color_hue() +
  geom_vline(xintercept = x)+
  theme_minimal()


example_2 =C[92, -c(1:3)] %>% as.numeric() %>% data.frame()
names(example_2) = "wTO"
example_2$run = 1:1000

delta = 0.2
min = (C[92,3] - delta) %>% as.numeric()
max = (C[92,3] + delta) %>% as.numeric()
example_2$p = ifelse( example_2$wTO > min & example_2$wTO < max, 0, 1)

example_2$prop = cumsum(example_2$p)
example_2$pval = example_2$prop/example_2$run


ggplot(example_2) +
  aes(x = run, y = pval) +
  geom_point(size = 1L, colour = "#bd3786") +
  geom_vline(xintercept = 100)+
  # geom_smooth(span = 1L, colour = "#bd3786") +
  theme_minimal()
