library(shrink)

# Load data ------------------------------------------------------
case1.bodyfat <-
  read.table("data/case1_bodyfat.txt", header = T, sep = ";")

# EPV --------------------------------------------------------------
pred <- c("age", "weight_kg", "height_cm", "neck", "chest", "abdomen", "hip", 
          "thigh", "knee", "ankle", "biceps", "forearm", "wrist")
epv <- dim(case1.bodyfat)[1] / length(pred)
epv

# Correlation structure --------------------------------------------
cor(case1.bodyfat[, pred])

# Estimate full model ----------------------------------------------
formula <- paste("siri~", paste(pred, collapse = "+"))
full_mod <- lm(formula, data = case1.bodyfat, x = T, y = T)
full_est <- coef(full_mod)
full_se <- coef(summary(full_mod))[, "Std. Error"]
summary(full_mod)

# Selected model ---------------------------------------------------
sel_mod <- step(lm(formula, data = case1.bodyfat,  x=T,y=T), 
                direction = "backward",
                scope = list(upper = formula, 
                             lower = formula(siri~abdomen+height_cm)),
                trace = 0)

summary(sel_mod)

sel_est <- coef(sel_mod)[c("(Intercept)", pred)]
sel_est[is.na(sel_est)] <- 0
names(sel_est) <- c("(Intercept)", pred)
sel_se <- coef(summary(sel_mod))[, "Std. Error"][c("(Intercept)", pred)]
sel_se[is.na(sel_se)] <- 0


# Bootstrap ----------------------------------------------------------
bootnum <- 1000
boot_est <-  boot_se <- matrix(0, ncol = length(pred) + 1, nrow = bootnum,
                               dimnames = list(NULL, c("(Intercept)", pred)))

set.seed(5437854)
for (i in 1:bootnum) {
  data_id <- sample(1:dim(case1.bodyfat)[1], replace = T)
  boot_mod <- step(lm(formula, data = case1.bodyfat[data_id, ], 
                             x=T, y=T), 
                         scope = list(upper = formula, 
                                      lower = 
                                        formula(siri ~ abdomen + height_cm)),
                         direction = "backward", trace = 0)
  boot_est[i, names(coef(boot_mod))] <- coef(boot_mod)
  boot_se[i, names(coef(boot_mod))] <- coef(summary(boot_mod))[, "Std. Error"]
}
boot_01 <- (boot_est != 0) * 1
boot_inclusion <- apply(boot_01, 2, function(x) sum(x) / length(x) * 100)



# Overview of estimates and measures --------------------------------
# Table 5 in the manuscript -----------------------------------------
sqe <- (t(boot_est) - full_est) ^ 2
rmsd <- apply(sqe, 1, function(x) sqrt(mean(x)))
rmsdratio <- rmsd / full_se
boot_mean <- apply(boot_est, 2, mean)
boot_meanratio <- boot_mean / full_est
boot_relbias <- (boot_meanratio / (boot_inclusion / 100) - 1) * 100
boot_median <- apply(boot_est, 2, median)
boot_025per <- apply(boot_est, 2, function(x) quantile(x, 0.025))
boot_975per <- apply(boot_est, 2, function(x) quantile(x, 0.975))
  
overview <- round(cbind(full_est, full_se, boot_inclusion, sel_est, sel_se, 
                        rmsdratio, boot_relbias, boot_median, boot_025per, 
                        boot_975per), 4)
overview <- overview[order(overview[,"boot_inclusion"], decreasing=T),]
overview


# Model frequency ---------------------------------------------------
# Table 6 in the manuscript -----------------------------------------
pred_ord <- names(boot_inclusion)[order(boot_inclusion, decreasing = T)]
boot_01 <- boot_01[, pred_ord]
boot_01 <- cbind(boot_01, count = rep(1, times = bootnum))
boot_modfreq <- aggregate(count ~ ., data = boot_01, sum)
boot_modfreq[, "percent"] <- boot_modfreq$count / bootnum * 100
boot_modfreq <- boot_modfreq[order(boot_modfreq[, "percent"], decreasing = T), ]
boot_modfreq[, "cum_percent"] <- cumsum(boot_modfreq$percent)
boot_modfreq <- boot_modfreq[boot_modfreq[, "cum_percent"] <= 80, ]
if (dim(boot_modfreq)[1] > 20) boot_modfreq <- boot_modfreq[1:20, ]


cbind("Predictors"= apply(boot_modfreq[,c(2:14)], 1, 
                          function(x) paste(names(x[x==1]), collapse=" ")),
      boot_modfreq[,c("count", "percent", "cum_percent")])


# Model frequency in % of selected model ----------------------------
sel_modfreq <- sum(apply(boot_01[, -dim(boot_01)[2]], 1, function(x)
    identical(((sel_est != 0) * 1)[pred_ord], x))) / bootnum * 100

sel_modfreq


# Pairwise inclusion frequency in % ----------------------------------
# Supplementary Table 2 in supporting information --------------------
pval <- 0.01
boot_pairfreq <- matrix(100, ncol = length(pred), nrow = length(pred),
                        dimnames = list(pred_ord[-1], 
                                        pred_ord[-1]))

expect_pairfreq <- NULL
combis <- combn(pred_ord[-1], 2)

for (i in 1:dim(combis)[2]) {
  boot_pairfreq[combis[1, i], combis[2, i]] <-
    sum(apply(boot_01[, combis[, i]], 1, sum) == 2) / bootnum * 100
  expect_pairfreq[i] <-
    boot_inclusion[grepl(combis[1, i], pred_ord)][1] *
    boot_inclusion[grepl(combis[2, i], pred_ord)][1] / 100
  boot_pairfreq[combis[2, i], combis[1, i]] <-
    ifelse(is(suppressWarnings(try(chisq.test(boot_01[, combis[1, i]], 
                                              boot_01[, combis[2, i]]), 
                                   silent = T)), "try-error"),
      NA, ifelse(suppressWarnings(
                   chisq.test(boot_01[, combis[1, i]],
                              boot_01[, combis[2, i]])$p.value) > pval,
        "", ifelse(as.numeric(boot_pairfreq[combis[1, i], combis[2, i]]) < 
                     expect_pairfreq[i], "-", "+")))
}

diag(boot_pairfreq) <- boot_inclusion[pred_ord][-1]

print(boot_pairfreq, quote = F)



# Shrinkage factors --------------------------------------------------
# Supplementary Table 3 in supporting information --------------------
sel_mod_shrinkg <- shrink(sel_mod, type="global")
sel_mod_shrinkp <- shrink(sel_mod, type="parameterwise")

sel_mod_shrinkg$ShrinkageFactors
sel_mod_shrinkp$ShrinkageFactors


sel_mod_shrinkp_vcov<-vcov(sel_mod_shrinkp)

round(cbind("Shrinkage factors" = sel_mod_shrinkp$ShrinkageFactors,
            "SE of shrinkage factors" = sqrt(diag(sel_mod_shrinkp_vcov))[-1],
            "Correlation matrix of shrinkage factors" = 
              cov2cor(sel_mod_shrinkp_vcov)[-1,-1]),4)


      