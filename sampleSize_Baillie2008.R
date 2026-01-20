## ====================== t2-only SRLI/GDI STATUS (LEVEL) ======================
## Goal: pick the smallest sample size n so the CURRENT indicator (at t2) is precise.
## This code does NOT use t1 and does NOT estimate trends.
##############################################
# RLI (t2-only) – full drop-in replacement
##############################################

## ---- Constants & helpers ----
RLI_WEIGHTS    <- c(LC=0, NT=1, VU=2, EN=3, CR=4, RE=5, EW=5, EX=5)
EXCLUDE_IN_RLI <- c("DD","NA")
RLI_MAXW       <- 5

norm_status <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("\\(.*\\)$", "", x)                 # e.g. "CR(PE)" -> "CR"
  ok <- c(names(RLI_WEIGHTS), EXCLUDE_IN_RLI)
  x[!(x %in% ok)] <- NA_character_
  x
}

pick_species_col <- function(df) {
  cand <- c("종명","종명(가나다순)","종명(국명)")
  hit  <- cand[cand %in% names(df)]
  if (!length(hit)) stop("No species-name column found: need one of ", paste(cand, collapse=", "))
  hit[1]
}

## Build t2 vector from your (already filtered) data.frame
## Keeps everything; DD/NA/NA are removed later inside rli() and bootstrap.
prepare_t2_vector <- function(df, species_col = NULL, t2_col = "개정판_평가범주") {
  if (is.null(species_col)) species_col <- pick_species_col(df)
  need <- c(species_col, t2_col)
  if (!all(need %in% names(df))) stop("Missing column(s): ", paste(setdiff(need, names(df)), collapse=", "))
  sp   <- trimws(as.character(df[[species_col]]))
  keep <- nzchar(sp)
  sp   <- make.names(sp[keep], unique = TRUE)
  t2   <- norm_status(df[[t2_col]][keep])
  names(t2) <- sp
  t2
}

## ---- RLI (excludes DD/NA and NA explicitly) ----
rli <- function(cats) {
  use <- !is.na(cats) & !(cats %in% EXCLUDE_IN_RLI)
  cats_use <- cats[use]
  if (!length(cats_use)) return(NA_real_)
  1 - sum(RLI_WEIGHTS[cats_use], na.rm = TRUE) / (RLI_MAXW * length(cats_use))
}

## ---- Bootstrap for t2-only stability (level test) ----
## Left axis metric = error rate: P(|RLI_subset - RLI_full| > tau_level)
bootstrap_rli_level <- function(cats_t2,
                                n_grid = NULL,
                                R = 20000,
                                tau_level = 0.05,   # tolerance band in RLI units
                                seed = 1,
                                step_n = 20,        # grid step (every 30 by default)
                                n_min  = 10,        # smallest n to try
                                ensure_nmax = TRUE) {
  set.seed(seed)
  ## eligible pool = non-missing and non-DD/NA
  pool <- which(!is.na(cats_t2) & !(cats_t2 %in% EXCLUDE_IN_RLI))
  Nmax <- length(pool)
  if (Nmax < 10L) stop("Too few species after excluding DD/NA (N < 10).")
  
  ## n-grid
  if (is.null(n_grid)) {
    if (Nmax < n_min) {
      start  <- max(5L, floor(Nmax/3))
      step   <- max(1L, floor(Nmax/10))
      n_grid <- seq(start, Nmax, by = step)
    } else {
      n_grid <- seq(n_min, Nmax, by = step_n)
    }
    if (ensure_nmax && (length(n_grid) == 0L || tail(n_grid, 1L) != Nmax)) {
      n_grid <- c(n_grid, Nmax)
    }
  }
  n_grid <- sort(unique(pmin(n_grid, Nmax)))
  n_grid <- n_grid[n_grid >= 2L]
  
  full_rli <- rli(cats_t2)
  
  out <- lapply(n_grid, function(n) {
    vals <- numeric(R)
    for (i in seq_len(R)) {
      idx    <- sample(pool, n, replace = FALSE)
      vals[i] <- rli(cats_t2[idx])
    }
    data.frame(
      n         = n,
      err_rate  = mean(abs(vals - full_rli) > tau_level),  # LEFT axis (probability)
      mean_rli  = mean(vals),
      mean_bias = mean(vals - full_rli),
      lo        = unname(quantile(vals, 0.025)),
      hi        = unname(quantile(vals, 0.975)),
      stringsAsFactors = FALSE
    )
  })
  bt <- do.call(rbind, out)
  attr(bt, "full_rli")  <- full_rli
  attr(bt, "tau_level") <- tau_level
  bt
}

## ---- Recommendation (t2-only) ----
## omega: target CI width; thr_pct: max allowed error rate; optional tau_bias gate
recommend_n_level <- function(bt,
                              omega = 0.05,
                              tau_level = NULL,
                              tau_bias = NA,
                              thr_pct = 0.05) {
  if (is.null(tau_level) && !is.null(attr(bt,"tau_level"))) tau_level <- attr(bt,"tau_level")
  bt$width <- bt$hi - bt$lo
  
  n_width <- bt$n[bt$width    <= omega][1]
  n_err   <- bt$n[bt$err_rate <= thr_pct][1]
  n_bias  <- if (!is.na(tau_bias)) bt$n[abs(bt$mean_bias) <= tau_bias][1] else NA_integer_
  
  n_rec <- max(n_width, n_err, n_bias, na.rm = TRUE)
  list(n_width = n_width, n_err = n_err, n_bias = n_bias, n_rec = n_rec)
}

## ---- Plot: left = error (%, optional), right = RLI (±95% CI) ----
plot_rli_level_two_axis <- function(bt,
                                    omega = 0.05,
                                    tau_level = NULL,
                                    n_rec = NA,
                                    thr_pct = 0.05,
                                    left_as_percent = TRUE,
                                    show_tau_band = FALSE,
                                    main = "RLI level stability (t2 only)") {
  full_rli <- attr(bt, "full_rli")
  if (is.null(tau_level) && !is.null(attr(bt,"tau_level"))) tau_level <- attr(bt,"tau_level")
  
  op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
  par(mar = c(4, 4, 2, 5))  # room for right axis label
  
  ## LEFT: error rate
  y_left <- bt$err_rate
  ylab_l <- "Error rate (|RLI_subset - RLI_full| > τ)"
  hline  <- thr_pct
  if (left_as_percent) { y_left <- 100*y_left; hline <- 100*thr_pct; ylab_l <- "Error rate (%)" }
  plot(bt$n, y_left, type="b", pch=16, xlab="Sample size (n)", ylab=ylab_l, main=main)
  abline(h = hline, lty = 2)
  
  ## RIGHT: RLI ±95% CI (with optional τ band)
  par(new=TRUE)
  yl <- range(c(bt$lo, bt$hi, full_rli), na.rm = TRUE)
  plot(bt$n, bt$mean_rli, type="n", axes=FALSE, xlab="", ylab="", ylim=yl)
  
  if (show_tau_band && !is.null(tau_level) && !is.na(tau_level)) {
    xr <- range(bt$n); dx <- 0.03*diff(xr)
    rect(xr[1]-dx, full_rli - tau_level, xr[2]+dx, full_rli + tau_level,
         col = adjustcolor("grey80", 0.45), border = NA)
  }
  
  width_ok <- (bt$hi - bt$lo) <= omega
  arrows(bt$n, bt$lo, bt$n, bt$hi, angle=90, code=3, length=0.04,
         col = ifelse(width_ok, "black", "grey65"))
  points(bt$n, bt$mean_rli, pch = ifelse(width_ok, 16, 1),
         col = ifelse(width_ok, "black", "grey35"))
  abline(h = full_rli, lty = 2, col = "grey30")
  axis(4); mtext("RLI (subset mean ± 95% CI)", side=4, line=3)
  
  if (!is.na(n_rec)) abline(v = n_rec, lty = 3, col = "red")
}

################# Example ######################
library(readxl)
# library(dplyr)  # if you want %>% filtering

raw <- read_excel("S:/2차년도_과제/보고서_논문_제안서_첨부/KR_RedList_1438_1530_shortList.xlsx", sheet = 1)

## Filter taxa (you control this step)
df_taxa  <- subset(raw, 분류군 %in% c("거미"))
df_taxa2<-df_taxa %>% filter(개정판_평가범주 %in% c("CR","RE","EN","VU","NT","LC"))
speciesNum<-nrow(df_taxa2) 
# TIP: t2-only 분석에는 '신규'도 유효한 t2 코드가 있으면 포함해도 됩니다.
# 만약 정말 빼고 싶다면:
# df_taxa <- subset(df_taxa, 변동현황 != "신규")

## Build the t2 vector
cats_t2 <- prepare_t2_vector(df_taxa2)  # or prepare_t2_vector(df_taxa, species_col="종명")

## Run bootstrap (every 30 by default) & recommend n
bt2  <- bootstrap_rli_level(cats_t2, R = 50000, tau_level = 0.05, seed = 1023, n_min  = 30,   step_n = 20)
rec  <- recommend_n_level(bt2, omega = NA, tau_level = 0.05, tau_bias = NA, thr_pct = 0.05)
rec;speciesNum

## Plot (two axes, with τ & ω guides)
plot_rli_level_two_axis(bt2, omega = 0.05, tau_level = NA, n_rec = rec$n_rec, thr_pct = 0.05)


140+139+230








