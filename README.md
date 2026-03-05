# ppmlen
elastic net for gravity analysis

# sample code using sample data
devtools::install_github("donginerickim/ppmlen")
library(ppmlen)


# load sample data
df <- read.csv(system.file("extdata","sample.csv",package="ppmlen"))
fit <- ppml_en_plugin(df, dep="trade_f", always_keep="tau")

# construct fixed effects
df$ik <- interaction(df$exp, df$itpd_id)
df$jk <- interaction(df$imp, df$itpd_id)
df$ij <- interaction(df$exp, df$imp)

# run PPML-EN
fit <- ppml_en_plugin(
  data = df,
  dep = "trade_f",
  fixed = list(c("ik","time"), c("jk","time"), c("ij","itpd_id")),
  cluster = c("ij"),
  alpha = 0.5,
  K = 2,
  always_keep = "tau"
)

fit$selected
