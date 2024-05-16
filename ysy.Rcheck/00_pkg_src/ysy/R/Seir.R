#' ODE solition
#' This function solves an ODE system based on the given parameter list.
#' @param user created parameters
#' @return solution of an ODE function system
#' @export
model <- function(t, y, param) {
  #设定参数
  S <- y[1]
  E <- y[2]
  I1 <- y[3]
  I2 <- y[4]
  R <- y[5]
  N <- param["N"]

  #设定参数
  beta1 <- param["beta1"]
  beta2 <- param["beta2"]
  mu1 <- param["mu1"]
  mu2 <- param["mu2"]
  lamda1 <- param["lamda1"]
  lamda2 <- param["lamda2"]
  gamma1 <- param["gamma1"]
  gamma2 <- param["gamma2"]
  mu <- param["mu"]

  #传染病数学模型（添加I1和I2）
  dSt <- mu * (N - S) - beta1 * S * I1/N - beta2 * S * I2/N
  dEt <- beta1 * S * I1/N + beta2 * S * I2/N - (mu + lamda1) * E
  dIt1 <- lamda1 * E - (mu1 + gamma1) * I1
  dIt2 <- lamda2 * E - (mu2 + gamma2) * I2
  dRt <- gamma1 * I1 + gamma2 * I2 - mu * R

  #求解结果整合成向量表达
  outcome <- c(dSt, dEt, dIt1, dIt2, dRt)

  #返回微分方程系统求解结果
  list(outcome)
}

#设置评估参数值
times <- seq(0, 156, by = 1/7)
param <- c(mu = 0.000, beta1 = 4, beta2 = 5, gamma1 = 0.01, gamma2 = 0.02, mu1 = 0.000, mu2 = 0.000, lamda1 = 0.004, lamda2 = 0.005, N = 1)
init <- c(S = 0.9999, E = 0.00008, I1 = 0.00001, I2 = 0.00005, R = 0)

#调用微分方程求解函数，传入初始条件，评估时间、模型以及参数信息
result <- deSolve::ode(y = init, times = times, func = model, parms = param)
result <- as.data.frame(result)

tail(round(result, 3.6), 10)

#结果画图
seirplot <- ggplot2::ggplot(data = result) +
  ggplot2::geom_line(ggplot2::aes(x = time, y = S, col = "S"), lwd = 2) +
  ggplot2::geom_line(ggplot2::aes(x = time, y = I1, col = "I1"), lwd = 2) +
  ggplot2::geom_line(ggplot2::aes(x = time, y = I2, col = "I2"), lwd = 2) +
  ggplot2::geom_line(ggplot2::aes(x = time, y = R, col = "R"), lwd = 2) +
  ggplot2::geom_line(ggplot2::aes(x = time, y = E, col = "E"), lwd = 2) +
  ggplot2::labs(x = "Time", y = "Ratio") +
  ggplot2::scale_color_manual(name = "SEIR",
                              values = c("S" = "orange", "E" = "purple", "I1" = "red", "I2" = "blue", "R" = "green"))

#绘制仿真结果并保存为矢量文件
seirplot
ggplot2::ggsave(seirplot, file = "seir.pdf", width = 7, height = 6)
ggplot2::ggsave(seirplot, file = "seir.svg", width = 7, height = 6)



