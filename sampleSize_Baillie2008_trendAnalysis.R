#### Code to calculate sample size from two time points(t1 and t2) using the Baillie et al. 2008 ####
# 1) previously I used the data (raw) that I extracted from the pdf.
# 2) I then, used the raw2 which I acquired from NIBR. Hence, the data names are different. 



## ====== Red List setup ======
RLI_WEIGHTS <- c(LC=0, NT=1, VU=2, EN=3, CR=4, RE=5, EW=5, EX=5)
EXCLUDE_IN_RLI <- c("DD","NA")   # excluded at each time in RLI
RLI_MAXW <- 5


## Extract English category code from Korean format like "최소관심(LC)" -> "LC"
extract_category <- function(x) {
  x <- as.character(x)
  # Extract text inside parentheses
  code <- sub(".*\\(([A-Z]+)\\).*", "\\1", x)
  # If no parentheses found, return NA
  code[!grepl("\\([A-Z]+\\)", x)] <- NA_character_
  code
}

## Pick species-name column present in your data.frame
pick_species_col <- function(df) {
  cand <- c("종명","종명(가나다순)","종명(국명)")
  hit <- cand[cand %in% names(df)]
  if (length(hit) == 0) stop("No species-name column found: need one of ", paste(cand, collapse=", "))
  hit[1]
}

## Convert a filtered data.frame -> named t1,t2 vectors (no filtering/IO here)
prepare_status_vectors <- function(df, species_col = NULL,
                                   t1_col = "초판_평가범주", t2_col = "개정판_평가범주") {
  if (is.null(species_col)) species_col <- pick_species_col(df)
  need <- c(species_col, t1_col, t2_col)
  if (!all(need %in% names(df))) stop("Missing columns: ", paste(setdiff(need, names(df)), collapse=", "))
  sp <- trimws(as.character(df[[species_col]]))
  keep <- nzchar(sp)
  sp <- make.names(sp[keep], unique = TRUE)
  t1 <- norm_status(df[[t1_col]][keep])
  t2 <- norm_status(df[[t2_col]][keep])
  names(t1) <- names(t2) <- sp
  list(t1 = t1, t2 = t2)
}

## ====== RLI & ΔRLI ======
rli <- function(cats) {
  use <- !(cats %in% EXCLUDE_IN_RLI)
  cats_use <- cats[use]
  if (!length(cats_use)) return(NA_real_)
  1 - sum(RLI_WEIGHTS[cats_use], na.rm = TRUE) / (RLI_MAXW * length(cats_use))
}

## ΔRLI between t1 and t2
## mode="overlap": only species non-DD/NA at BOTH times (closest to Baillie)
## mode="timepoint-specific": exclude DD/NA separately at each time
delta_rli <- function(cats_t1, cats_t2, mode = c("overlap","timepoint-specific")) {
  mode <- match.arg(mode)
  if (mode == "overlap") {
    keep <- !(cats_t1 %in% EXCLUDE_IN_RLI) & !(cats_t2 %in% EXCLUDE_IN_RLI)
    r1 <- rli(cats_t1[keep]); r2 <- rli(cats_t2[keep])
  } else {
    r1 <- rli(cats_t1);        r2 <- rli(cats_t2)
  }
  list(delta = unname(r2 - r1),
       n_used_t1 = sum(!(cats_t1 %in% EXCLUDE_IN_RLI)),
       n_used_t2 = sum(!(cats_t2 %in% EXCLUDE_IN_RLI)),
       n_overlap = sum(!(cats_t1 %in% EXCLUDE_IN_RLI) & !(cats_t2 %in% EXCLUDE_IN_RLI)))
}

## ====== Bootstrap SRLI across sample sizes (Baillie sign test) ======
bootstrap_srli <- function(cats_t1, cats_t2,
                           n_grid = NULL,
                           R = 20000,
                           mode = "overlap",
                           seed = 1,
                           step_n = 30,     # <- step size between sample sizes
                           n_min  = 30,     # <- smallest n to try (default 30)
                           ensure_nmax = TRUE) {
  set.seed(seed)
  
  pool <- if (mode == "overlap") {
    which(!(cats_t1 %in% EXCLUDE_IN_RLI) & !(cats_t2 %in% EXCLUDE_IN_RLI))
  } else {
    seq_along(cats_t1)
  }
  Nmax <- length(pool)
  if (Nmax < 10L) stop("Too few species in pool to bootstrap (N<10).")
  
  ## ----- make the n grid -----
  if (is.null(n_grid)) {
    if (Nmax < n_min) {
      # fallback for very small pools: try a small grid from ~Nmax/3 to Nmax
      start <- max(10L, floor(Nmax/3))
      step  <- max(5L, floor(Nmax/10))
      n_grid <- seq(start, Nmax, by = step)
    } else {
      n_grid <- seq(n_min, Nmax, by = step_n)
    }
    if (ensure_nmax && (length(n_grid) == 0L || tail(n_grid, 1L) != Nmax)) {
      n_grid <- c(n_grid, Nmax)
    }
  }
  n_grid <- sort(unique(pmin(n_grid, Nmax)))
  n_grid <- n_grid[n_grid >= 2L]  # sanity
  
  full_delta <- delta_rli(cats_t1, cats_t2, mode)$delta
  
  out <- lapply(n_grid, function(n) {
    deltas <- numeric(R)
    for (i in seq_len(R)) {
      idx <- sample(pool, n, replace = FALSE)
      deltas[i] <- delta_rli(cats_t1[idx], cats_t2[idx], mode = "overlap")$delta
    }
    data.frame(
      n          = n,
      wrong_dir  = mean(sign(deltas) != sign(full_delta)),
      mean_delta = mean(deltas),
      mean_bias  = mean(deltas - full_delta),
      lo         = unname(quantile(deltas, 0.025)),
      hi         = unname(quantile(deltas, 0.975)),
      stringsAsFactors = FALSE
    )
  })
  bt <- do.call(rbind, out)
  attr(bt, "full_delta") <- full_delta
  bt
}


## ====== Plotters (Fig-1 style) ======
plot_srli_single_panel <- function(bt,
                                   wrong_thresh = 0.05,
                                   left_as_percent = TRUE,
                                   main = "SRLI trend test: sign error & ΔRLI") {
  full_delta <- attr(bt, "full_delta")
  
  op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
  par(mar = c(4, 4, 2, 5))  # <- extra right margin so the label is visible
  
  ## LEFT axis: wrong-direction (probability or %)
  y_left <- bt$wrong_dir
  ylab_left <- "Wrong-direction probability"
  hline <- wrong_thresh
  if (left_as_percent) {
    y_left <- 100 * y_left
    ylab_left <- "Wrong-direction (%)"
    hline <- 100 * wrong_thresh
  }
  plot(bt$n, y_left, type = "b", pch = 16,
       xlab = "Sample size (n)", ylab = ylab_left, main = main)
  abline(h = hline, lty = 2)
  
  ## RIGHT axis: ΔRLI (subset mean ± 95% CI), horizontal line = FULL ΔRLI
  yl <- range(c(bt$lo, bt$hi, full_delta), na.rm = TRUE)
  par(new = TRUE)
  plot(bt$n, bt$mean_delta, type = "n", axes = FALSE, xlab = "", ylab = "", ylim = yl)
  arrows(bt$n, bt$lo, bt$n, bt$hi, angle = 90, code = 3, length = 0.04, col = "grey40")
  points(bt$n, bt$mean_delta, pch = 1, col = "grey10")
  abline(h = full_delta, lty = 2, col = "grey30")
  axis(4)
  mtext("ΔRLI (mean ± 95% CI, t2 − t1)", side = 4, line = 3)  # <- explicit right label
}

plot_srli_two_panel <- function(bt, wrong_thresh = 0.05) {
  full_delta <- attr(bt, "full_delta")
  op <- par(mfrow=c(1,2), mar=c(4,4,2,1)); on.exit(par(op), add=TRUE)
  plot(bt$n, bt$wrong_dir, type="b", pch=16,
       xlab="Sample size (n)", ylab="Wrong-direction probability",
       main="Trend sign error vs sample size")
  abline(h = wrong_thresh, lty = 2)
  
  yl <- range(c(bt$lo, bt$hi, full_delta), na.rm = TRUE)
  plot(bt$n, bt$mean_delta, type="b", pch=16, ylim=yl,
       xlab="Sample size (n)", ylab="ΔRLI (mean ± 95% CI)",
       main="ΔRLI vs sample size")
  arrows(bt$n, bt$lo, bt$n, bt$hi, angle=90, code=3, length=0.04)
  abline(h = full_delta, lty = 2)
}

recommend_n_sign <- function(bt, wrong_thresh = 0.05) {
  idx <- which(bt$wrong_dir <= wrong_thresh)
  if (length(idx)) bt$n[min(idx)] else NA_integer_
}


################## Example ####################
library(dplyr)
library(readxl)
## You do the I/O
raw <- read_excel("S:/2차년도_과제/보고서_논문_제안서_첨부/KR_RedList_1438_1530_shortList.xlsx", sheet = 1)
raw2 <- read_excel("S:/2차년도_과제/보고서_논문_제안서_첨부/국가생물적색목록_평가현황(24.12월_기준).xlsx", sheet = 1)
str(raw2)

## Clean up the data for the ones that have both t1 and t2; 
str(raw2.1)
raw2.1<- raw2 %>% select(no., 분류군, 국명,평가년도...5, 평가범주...6, 평가년도...8, 평가범주...9, 변동현황,멸종위기종) %>%
         mutate(초판_년도=평가년도...5, 
                초판_평가범주=평가범주...6,
                개정판_년도=평가년도...8, 
                개정판_평가범주=평가범주...9,
                종명=국명) %>%
         mutate(초판_평가범주 = extract_category(초판_평가범주),
          개정판_평가범주 = extract_category(개정판_평가범주),
          초판_년도=as.numeric(초판_년도), 개정판_년도= as.numeric(개정판_년도)
         )  %>% 
  select(no., 분류군, 종명, 초판_년도,초판_평가범주,개정판_년도, 개정판_평가범주, 변동현황) %>%
        filter(초판_년도 > 2000 & 개정판_년도 > 2000)
        
  

  

## View sample size
categoryGivenSamples<-raw2.1 %>% select(분류군)  %>% table()
categoryGivenSamplesNotDD<-raw2.1 %>% filter(변동현황 != "신규" & !is.na(변동현황) & 
                    !is.na(초판_평가범주) & !is.na(개정판_평가범주) &
                    초판_평가범주 !="DD" & 개정판_평가범주 !="DD" ) %>%
                    select(분류군)  %>% table()


## You do the filtering (example: plants & insects)
df_taxa <- subset(raw2.1, 분류군 %in% c("관속식물"))  # change taxa list as needed
df_taxa
df_taxa2<-df_taxa %>% filter(변동현황 != "신규" & !is.na(변동현황) & 
                               !is.na(초판_평가범주) & !is.na(개정판_평가범주) &
                            초판_평가범주 !="DD" & 개정판_평가범주 !="DD" ) 
df_taxa2
# df_taxa2 %>% filter(분류군=="조류") %>% data.frame()

## Build t1/t2 vectors from your filtered df
v  <- prepare_status_vectors(df_taxa2)  # or specify species_col=... if needed
t1 <- v$t1
t2 <- v$t2

## Compute the full (truth) ΔRLI for the interval
truth <- delta_rli(t1, t2, mode = "overlap")
truth$delta; truth$n_overlap

## Bootstrap across sample sizes and plot
bt <- bootstrap_srli(t1, t2, R = 50000, mode = "overlap",  n_min  = 20, step_n = 10,seed = 1023)
plot_srli_single_panel(bt, wrong_thresh = 0.05)   # Fig-1 style
# or:
# plot_srli_two_panel(bt, wrong_thresh = 0.05)

## Recommended sample size by Baillie 5% sign rule
n_rec <- recommend_n_sign(bt, wrong_thresh = 0.05)
n_rec


# 국가적색목록 중 평가범주 부여된 종
95+39+64+34+416+1822+389+339+290+689
# 국가적색목록 중 초판 개정판 모두 있는 종
95+36+60+31+372+737+167+289+210+564
#추세분석 결과 최소 종 수
40+35+50+20+30+460+120+160+120+40
