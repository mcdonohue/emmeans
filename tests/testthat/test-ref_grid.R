context("Reference grids")

pigs.lm = lm(log(conc) ~ source + factor(percent), data = pigs)
rg = ref_grid(pigs.lm)
rg1 = ref_grid(pigs.lm, at = list(source = "soy", percent = 12))

pigs.lm2 = update(pigs.lm, conc ~ source + percent, data = pigs)
rg2 = ref_grid(pigs.lm2)
rg2a = ref_grid(pigs.lm2, at = list(source = c("fish", "soy"), percent = 10))
rg2c = ref_grid(pigs.lm2, cov.reduce = FALSE)
rg2m = ref_grid(pigs.lm2, cov.reduce = min)

pigs.lm3 = update(pigs.lm2, . ~ source + source:factor(percent))
pigs = transform(pigs, sp = interaction(source, percent))
pigs.lm4 = update(pigs.lm2, . ~ source + sp)

test_that("Reference grid is constructed correctly", {
    expect_equal(nrow(rg@grid), 12)
    expect_equal(nrow(rg1@grid), 1)
    expect_equal(nrow(rg2@grid), 3)
    expect_equal(nrow(rg2a@grid), 2)
    expect_equal(nrow(rg2c@grid), 12)
    expect_equal(nrow(rg1@grid), 1)
    expect_equal(length(rg@levels), 2)
    expect_equal(rg2@levels$percent, mean(pigs$percent))
    expect_equal(rg2m@levels$percent, min(pigs$percent))
})

test_that("Reference grid extras are detected", {
    expect_equal(rg@misc$tran, "log")
    expect_true(is.null(rg2@misc$tran))
    expect_true(is.null(rg2@model.info$nesting))
    expect_is(ref_grid(pigs.lm3)@model.info$nesting, "list") # see note above
    expect_is(ref_grid(pigs.lm4)@model.info$nesting, "list") # see note above
})

colnames(ToothGrowth) <- c('len', 'choice of supplement', 'dose')
model <- stats::aov(`len` ~ `choice of supplement`, ToothGrowth)

test_that("Reference grid handles variables with spaces", {
    expect_output(str(ref_grid(model, ~`choice of supplement`)), "choice of supplement")
})

# models outside of data.frames
x = 1:10
y = rnorm(10)
mod1 = with(pigs, lm(log(conc) ~ source + factor(percent)))
test_that("ref_grid works with no data or subset", {
    expect_silent(ref_grid(lm(y ~ x)))
    expect_silent(ref_grid(mod1))
})

context("Tests for gls, lme, lmer")
library(pbkrtest)
library(nlme)
library(lme4)

fm1.lmer <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
fm1.lme <- lme(Reaction ~ Days, sleepstudy, random=~1 | Subject)
fm1.gls <- gls(Reaction ~ Days, sleepstudy, correlation = corCompSymm(form = ~ 1 | Subject))

rg.lmer = ref_grid(fm1.lmer)
rg.lme = ref_grid(fm1.lme)
rg.gls = ref_grid(fm1.gls)

test_that("Reference grid is constructed correctly for gls, lme, lmer.", {
    expect_equal(nrow(rg.lmer@grid), nrow(rg.lme@grid))
    expect_equal(nrow(rg.lmer@grid), nrow(rg.gls@grid))
    expect_equal(length(rg.lmer@levels), length(rg.lme@levels))
    expect_equal(length(rg.lmer@levels), length(rg.gls@levels))
    expect_equal(rg.lmer@levels$Days, rg.lme@levels$Days)
    expect_equal(rg.lmer@levels$Days, rg.gls@levels$Days)
    expect_equal(
        attr(as.data.frame(emmeans(fm1.lmer, specs=c('Days'))), 'mesg'),
        attr(as.data.frame(emmeans(fm1.lme, specs=c('Days'))), 'mesg'))
    expect_equal(
        attr(as.data.frame(emmeans(fm1.lmer, specs=c('Days'))), 'mesg'),
        attr(as.data.frame(emmeans(fm1.gls, specs=c('Days'))), 'mesg'))
    expect_equal(as.data.frame(emmeans(fm1.lmer, specs=c('Days'))), 
        as.data.frame(emmeans(fm1.lme, specs=c('Days'))))
    expect_equal(as.data.frame(emmeans(fm1.lmer, specs=c('Days'))), 
        as.data.frame(emmeans(fm1.gls, specs=c('Days'))), tolerance=0.00001)
    expect_equal(as.data.frame(emmeans(fm1.lmer, specs=c('Days')))$df, 
        as.data.frame(emmeans(fm1.gls, specs=c('Days')))$df)
})


data(Orthodont)
Orth1 <- Orthodont
Orth1$visno <- rep(1:4, nrow(Orth1)/4)
set.seed(20181011)
Orth1 <- Orth1[sample(0:1, size=nrow(Orth1), replace=TRUE, prob=c(0.75, 0.25)) == 0, ]
Orth2 <- Orth1[sample(1:nrow(Orth1)), ]
# so getData.gls can find it:
assign('Orth1', Orth1, envir=.GlobalEnv)
assign('Orth2', Orth2, envir=.GlobalEnv)

fm01.lmer <- lmer(distance ~ 1 + (age|Subject), data = Orth1)
fm11.lmer <- lmer(distance ~ age + Sex + (age|Subject), data = Orth1)
fm02.lmer <- lmer(distance ~ 1 + (age|Subject), data = Orth2)
fm12.lmer <- lmer(distance ~ age + Sex + (age|Subject), data = Orth2)
kr1.lmer <- KRmodcomp(fm11.lmer, fm01.lmer)
kr2.lmer <- KRmodcomp(fm12.lmer, fm02.lmer)
expect_equal(kr1.lmer$test, kr2.lmer$test, tolerance=1e-05)

fm01.lme <- lme(distance ~ 1, data = Orth1, random = ~ age|Subject)
fm11.lme <- lme(distance ~ age + Sex, data = Orth1, random = ~ age|Subject)
fm02.lme <- lme(distance ~ 1, data = Orth2, random = ~ age|Subject)
fm12.lme <- lme(distance ~ age + Sex, data = Orth2, random = ~ age|Subject)
kr1.lme <- KRmodcomp(fm11.lme, fm01.lme)
kr2.lme <- KRmodcomp(fm12.lme, fm02.lme)
expect_equal(kr1.lme$test, kr2.lme$test, tolerance=1e-05)

