# 测试 R/Component.R —— 文档见 docs/component.typ
# 四个函数统一用 ::: 取用，与 NAMESPACE 的导出状态解耦。
RemoveSeasonal <- rpkg.template:::RemoveSeasonal
lm_solve       <- rpkg.template:::lm_solve
LS_OneBreak    <- rpkg.template:::LS_OneBreak
LS_MultiBreaks <- rpkg.template:::LS_MultiBreaks

# RemoveSeasonal ---------------------------------------------------------------

test_that("RemoveSeasonal: 分月日均值正确", {
  y <- c(1, 3, 2, 6)
  t_mmdd <- c(101L, 101L, 102L, 102L)

  out <- RemoveSeasonal(y, t_mmdd)

  expect_named(out, c("mu", "anomaly"))
  expect_equal(out$mu, c(2, 4))
  expect_equal(out$anomaly, c(-1, 1, -2, 2))
})

test_that("RemoveSeasonal: mu 按 MMDD 升序对齐，与输入顺序无关", {
  y <- c(1, 3, 2, 6)
  t_mmdd <- c(102L, 101L, 102L, 101L)

  out <- RemoveSeasonal(y, t_mmdd)

  # 101 -> mean(3, 6) = 4.5；102 -> mean(1, 2) = 1.5
  expect_equal(out$mu, c(4.5, 1.5))
  expect_equal(out$anomaly, c(-0.5, -1.5, 0.5, 1.5))
})

test_that("RemoveSeasonal: 各组异常值均值为 0，长度与 y 一致", {
  set.seed(42)
  t_mmdd <- rep(seq_len(366L), times = 3L)
  y <- rnorm(length(t_mmdd)) + 5 * sin(2 * pi * t_mmdd / 366)

  out <- RemoveSeasonal(y, t_mmdd)

  expect_length(out$anomaly, length(y))
  expect_length(out$mu, 366L)
  expect_equal(
    tapply(out$anomaly, t_mmdd, mean),
    setNames(numeric(366L), seq_len(366L))
  )
})

test_that("RemoveSeasonal: y 中的 NA 不参与均值", {
  out <- RemoveSeasonal(c(1, NA, 3), c(101L, 101L, 101L))

  expect_equal(out$mu, 2)
  expect_equal(out$anomaly, c(-1, NA, 1))
})

test_that("RemoveSeasonal: 单元素组退化为异常值 0", {
  out <- RemoveSeasonal(c(5, 7), c(101L, 102L))

  expect_equal(out$mu, c(5, 7))
  expect_equal(out$anomaly, c(0, 0))
})

# B2：t_mmdd 为 NA 的位置不进入任何分组，当前实现静默返回 0。
# 此处只锁定"不报错、长度不变"，不锁定 0 这一错误行为，待修复后再加断言。
test_that("RemoveSeasonal: t_mmdd 含 NA 不报错且保持长度", {
  out <- RemoveSeasonal(c(1, 2, 3), c(101L, NA, 101L))

  expect_length(out$anomaly, 3L)
})

# lm_solve ---------------------------------------------------------------------

test_that("lm_solve: 拟合值加残差等于响应值，SSE 等于残差平方和", {
  set.seed(1)
  X <- seq_len(20)
  Y <- 1 + 0.5 * X + rnorm(20, sd = 0.01)

  out <- lm_solve(cbind(1, X), Y)

  expect_equal(out$fitted + out$resi, Y)
  expect_equal(out$SSE, sum(out$resi^2))
})

test_that("lm_solve: only.SSE 只返回 SSE", {
  X <- seq_len(10)
  Y <- 2 * X

  out <- lm_solve(cbind(1, X), Y, only.SSE = TRUE)

  expect_named(out, "SSE")
  expect_lt(out$SSE, 1e-20)
})

# LS_OneBreak ------------------------------------------------------------------

test_that("LS_OneBreak: 无噪声两段数据精确还原截距跳变", {
  X <- seq_len(20)
  Y <- ifelse(X <= 10, 1 + 0.5 * X, 3 + 0.5 * X)

  out <- LS_OneBreak(Y, X, 10L)

  expect_equal(unname(out$coef), c(1, 0.5, 2), tolerance = 1e-10)
  expect_equal(out$fitted, Y, tolerance = 1e-10)
  expect_lt(out$SSE, 1e-20)
})

test_that("LS_OneBreak: 断点模型 SSE 不大于无断点模型", {
  set.seed(1)
  X <- seq_len(60)
  Y <- ifelse(X <= 30, 0.2 * X, 5 + 0.2 * X) + rnorm(60, sd = 0.1)

  expect_lt(LS_OneBreak(Y, X, 30L)$SSE, lm_solve(cbind(1, X), Y)$SSE)
})

test_that("LS_OneBreak: Ic = NA 退化为线性回归", {
  X <- seq_len(10)
  Y <- 1 + 2 * X

  out <- LS_OneBreak(Y, X, NA)

  expect_length(out$coef, 2L)
  expect_equal(unname(out$coef), c(1, 2), tolerance = 1e-10)
})

# B4：断点落在边界时示性列与截距列共线，秩亏只以 NA 系数体现。
test_that("LS_OneBreak: 断点取边界时秩亏，系数含 NA", {
  X <- seq_len(10)
  Y <- 1 + 0.5 * X

  expect_true(anyNA(LS_OneBreak(Y, X, 0L)$coef))  # 示性列全为 1
  expect_true(anyNA(LS_OneBreak(Y, X, 10L)$coef)) # 示性列全为 0
})

# LS_MultiBreaks ---------------------------------------------------------------

test_that("LS_MultiBreaks: 三段数据还原共用斜率与截距跳变", {
  X <- seq_len(30)
  Y <- c(0 + 0.2 * X[1:10], 1 + 0.2 * X[11:20], 3 + 0.2 * X[21:30])

  out <- LS_MultiBreaks(Y, X, c(10L, 20L, 30L))

  # X 已中心化：beta0 为段 1 在 X = mean(X) 处的水平
  expect_equal(unname(out$coef), c(0.2 * mean(X), 0.2, 1, 3), tolerance = 1e-10)
  expect_equal(out$fitted, Y, tolerance = 1e-10)
})

test_that("LS_MultiBreaks: 单断点与 LS_OneBreak 拟合值、SSE 一致", {
  set.seed(2)
  X <- seq_len(40)
  Y <- ifelse(X <= 15, 1 + 0.3 * X, 2.5 + 0.3 * X) + rnorm(40, sd = 0.1)

  one <- LS_OneBreak(Y, X, 15L)
  multi <- LS_MultiBreaks(Y, X, c(15L, 40L))

  expect_equal(one$fitted, multi$fitted, tolerance = 1e-10)
  expect_equal(one$SSE, multi$SSE, tolerance = 1e-10)
})

test_that("LS_MultiBreaks: Ips 只有末位时退化为线性回归", {
  X <- seq_len(12)
  Y <- 3 - 0.5 * X

  out <- LS_MultiBreaks(Y, X, 12L)

  expect_length(out$coef, 2L)
  expect_equal(unname(out$coef), c(3 - 0.5 * mean(X), -0.5), tolerance = 1e-10)
})
