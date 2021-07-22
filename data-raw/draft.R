rm(list = ls())
library(matlab)
library(rootSolve)
library(tidyverse)
library(ggplot2)
library(patchwork)
devtools::load_all()
source('E:/17_HuaDong/learn/HM2009/RCode/R/ramsey3_1.R')
# para setting
alpha <- 0.27
betax <- 0.994
tmax <- 60

x0 <- ones(tmax,1)
bounds <- data.frame(dw =ones(tmax,1)*0.01, up = ones(tmax,1)*30)

# draw
draw_root <- function(ans){
  picdata <- data.frame(tt = 1:length(ans$root), k = ans$root)
  return(ggplot(picdata, aes(x = tt, y = k)) + geom_line())
}

# 非线性方程组求解
# 解这样得方程组，初值得选择真是太重要了
# dfsane(par = as.numeric(x0)-0.3, fn = ramsey,k0 = 0.02, kT = 0, delta = 1)

FixVMN2(x0 = x0,bounds = bounds, f = ramsey, k0 = 4.4, kT = 0, delta = 0.011)

ans <- multiroot(f = ramsey,start = x0*0.2-0.1,k0 = 0.02, kT = 0, delta = 1)
p1 <- draw_root(ans)
ans <- multiroot(f = ramsey,start = x0,k0 = 66.06, kT = 0, delta = 0.011)
p2 <- draw_root(ans)
ans <- multiroot(f = ramsey,start = x0-0.5,k0 = 0.25, kT = 0, delta = 1)
p3 <- draw_root(ans)
ans <- multiroot(f = ramsey,start = x0,k0 = 4.4, kT = 0, delta = 0.011)
p4 <- draw_root(ans)
(p1+p4)/(p3+p2)


# 观察初值是否会使得消费小于0
# ramsey(k = x0,k0 = 66.06,kT =0, eta = 1,delta = 0.011)
