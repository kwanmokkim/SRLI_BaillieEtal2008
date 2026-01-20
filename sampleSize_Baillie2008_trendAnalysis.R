#### Code to calculate sample size from two time points(t1 and t2) using the Baillie et al. 2008 ####
# 1) previously I used the data (raw) that I extracted from the pdf.
# 2) I then, used the raw2 which I acquired from NIBR. Hence, the data names are different.
# 3) If you have a tidy data run the functions and then move to step #2 for analyses. 


# This code implements the Sampled Red List Index (SRLI) methodology from:
# Baillie, J.E.M., et al. (2008). Toward monitoring global biodiversity. 
# Conservation Letters, 1(1), 18-26.


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
  cand <- c("종명","종명(가나다순)","종명(국명)","species","speciesName")
  hit <- cand[cand %in% names(df)]
  if (length(hit) == 0) stop("No species-name column found: need one of ", paste(cand, collapse=", "))
  hit[1]
}

norm_status <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("\\(.*\\)$", "", x)                 # e.g. "CR(PE)" -> "CR"
  ok <- c(names(RLI_WEIGHTS), EXCLUDE_IN_RLI)
  x[!(x %in% ok)] <- NA_character_
  x
}

## Convert a filtered data.frame -> named t1,t2 vectors (no filtering/IO here)
prepare_status_vectors <- function(df, species_col = NULL,
                                   t1_col = "early_category", t2_col = "later_category") {
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


## ====== Find exact n for threshold using binary search ======
## This function finds the precise sample size where error rate crosses the threshold,
## rather than being limited by step_n resolution.
##
## Parameters:
##   cats_t1, cats_t2: named vectors of IUCN categories at time 1 and 2
##   wrong_thresh: target threshold (default 0.05 = 5% as per Baillie et al.)
##   R: number of bootstrap replicates (default 50000)
##   mode: "overlap" (recommended) or "timepoint-specific"
##   seed: random seed for reproducibility
##   tol: stop binary search when range <= tol (default 1 = exact integer)
##   verbose: print progress messages
##
## Returns a list with:
##   n_exact: the minimum sample size where error rate <= threshold
##   final_rate: the actual error rate at n_exact
##   threshold: the threshold used
##   full_delta: the true ΔRLI for the full dataset
##   Nmax: maximum possible sample size

find_exact_n_threshold <- function(cats_t1, cats_t2,
                                   wrong_thresh = 0.05,
                                   R = 50000,
                                   mode = "overlap",
                                   seed = 1,
                                   tol = 1,
                                   verbose = TRUE) {
  set.seed(seed)
  
  # Get the pool of valid species
  pool <- if (mode == "overlap") {
    which(!(cats_t1 %in% EXCLUDE_IN_RLI) & !(cats_t2 %in% EXCLUDE_IN_RLI))
  } else {
    seq_along(cats_t1)
  }
  Nmax <- length(pool)
  
  if (Nmax < 10L) stop("Too few species in pool to search (N<10).")
  
  
  # Full dataset delta (truth)
  full_delta <- delta_rli(cats_t1, cats_t2, mode)$delta
  
  # Helper function: compute wrong-direction rate for a given n
  compute_wrong_rate <- function(n, seed_offset = 0) {
    set.seed(seed + seed_offset)
    deltas <- numeric(R)
    for (i in seq_len(R)) {
      idx <- sample(pool, n, replace = FALSE)
      deltas[i] <- delta_rli(cats_t1[idx], cats_t2[idx], mode = "overlap")$delta
    }
    mean(sign(deltas) != sign(full_delta))
  }
  
  # Step 1: Coarse search to find initial bounds
  if (verbose) cat("Step 1: Coarse search to find bounds...\n")
  
  # Create a coarse grid (about 20 points across the range)
  coarse_step <- max(10, floor(Nmax / 20))
  coarse_grid <- seq(20, Nmax, by = coarse_step)
  if (tail(coarse_grid, 1) != Nmax) coarse_grid <- c(coarse_grid, Nmax)
  
  # Ensure we start small enough
  if (coarse_grid[1] > 20) coarse_grid <- c(20, coarse_grid)
  coarse_grid <- coarse_grid[coarse_grid <= Nmax]
  
  lower <- min(coarse_grid)
  upper <- Nmax
  found_upper <- FALSE
  
  for (i in seq_along(coarse_grid)) {
    n <- coarse_grid[i]
    rate <- compute_wrong_rate(n, seed_offset = i)
    if (verbose) cat(sprintf("  n = %d, wrong_rate = %.4f\n", n, rate))
    
    if (rate <= wrong_thresh) {
      upper <- n
      found_upper <- TRUE
      # Set lower to previous grid point (or minimum)
      if (i > 1) lower <- coarse_grid[i - 1]
      break
    } else {
      lower <- n
    }
  }
  
  # If even at Nmax we exceed threshold, return Nmax with warning
  if (!found_upper) {
    final_rate <- compute_wrong_rate(Nmax, seed_offset = 999)
    if (final_rate > wrong_thresh) {
      warning("Even at maximum sample size (", Nmax, "), error rate (", 
              round(final_rate, 4), ") exceeds threshold (", wrong_thresh, ").")
      return(list(
        n_exact = Nmax,
        final_rate = final_rate,
        threshold = wrong_thresh,
        full_delta = full_delta,
        Nmax = Nmax,
        converged = FALSE
      ))
    } else {
      upper <- Nmax
    }
  }
  
  # Step 2: Binary search for exact threshold
  if (verbose) cat(sprintf("\nStep 2: Binary search between %d and %d...\n", lower, upper))
  
  iteration <- 0
  while ((upper - lower) > tol) {
    iteration <- iteration + 1
    mid <- floor((lower + upper) / 2)
    rate <- compute_wrong_rate(mid, seed_offset = 1000 + iteration)
    if (verbose) cat(sprintf("  n = %d, wrong_rate = %.4f\n", mid, rate))
    
    if (rate <= wrong_thresh) {
      upper <- mid
    } else {
      lower <- mid
    }
  }
  
  # Final verification at upper bound (with fresh seed)
  final_rate <- compute_wrong_rate(upper, seed_offset = 9999)
  
  # Due to stochastic variation, verify the result
  # If final_rate > threshold, increment by 1 and check again
  verify_attempts <- 0
  while (final_rate > wrong_thresh && upper < Nmax && verify_attempts < 5) {
    verify_attempts <- verify_attempts + 1
    upper <- upper + 1
    final_rate <- compute_wrong_rate(upper, seed_offset = 9999 + verify_attempts)
    if (verbose) cat(sprintf("  Verification: n = %d, wrong_rate = %.4f\n", upper, final_rate))
  }
  
  if (verbose) cat(sprintf("\n==> Result: n = %d (wrong_rate = %.4f)\n", upper, final_rate))
  
  list(
    n_exact = upper,
    final_rate = final_rate,
    threshold = wrong_thresh,
    full_delta = full_delta,
    Nmax = Nmax,
    converged = TRUE
  )
}


## ====== Plotters (Fig-1 style) ======
plot_srli_single_panel <- function(bt,
                                   wrong_thresh = 0.05,
                                   left_as_percent = TRUE,
                                   main = "SRLI trend test: sign error & ΔRLI",
                                   xlim = NULL,
                                   ylim_left = NULL,
                                   ylim_right = NULL) {
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
  
  ## Set default axis ranges if not provided
  if (is.null(xlim)) xlim <- range(bt$n)
  if (is.null(ylim_left)) ylim_left <- c(0, max(y_left, na.rm = TRUE) * 1.1)
  if (is.null(ylim_right)) ylim_right <- range(c(bt$lo, bt$hi, full_delta), na.rm = TRUE)
  
  ## Plot left axis (wrong-direction)
  plot(bt$n, y_left, type = "b", pch = 16,
       xlab = "Sample size (n)", ylab = ylab_left, main = main,
       xlim = xlim, ylim = ylim_left)
  abline(h = hline, lty = 2)
  
  ## RIGHT axis: ΔRLI (subset mean ± 95% CI), horizontal line = FULL ΔRLI
  par(new = TRUE)
  plot(bt$n, bt$mean_delta, type = "n", axes = FALSE, xlab = "", ylab = "", 
       xlim = xlim, ylim = ylim_right)
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


## ====== Recommend n from grid (original method) ======
## Returns the first n in the bootstrap grid where wrong_dir <= threshold
## Note: Resolution is limited by step_n used in bootstrap_srli()
## For exact threshold, use find_exact_n_threshold() instead

recommend_n_sign <- function(bt, wrong_thresh = 0.05) {
  idx <- which(bt$wrong_dir <= wrong_thresh)
  if (length(idx)) bt$n[min(idx)] else NA_integer_
}


################## Example ####################
# modify the code below to run the analysis
# If your code is already tidy go to step 2 for analyses. 

 library(dplyr)
 library(readxl)
# 
# #### 1. Import & Clean data ######
# # The codes are relevant for South Korea data. It is better to clean your data yourself.
# 
# ## You do the I/O
raw <- read_excel("S:/2차년도_과제/보고서_논문_제안서_첨부/KR_RedList_1438_1530_shortList.xlsx", sheet = 1)
raw2 <- read_excel("S:/2차년도_과제/보고서_논문_제안서_첨부/국가생물적색목록_평가현황(24.12월_기준).xlsx", sheet = 1)
str(raw2)

 
# ## Clean up the data; with both t1 and t2
# ## Your data should have columns for:
# ##   - species name
# ##   - early_category (IUCN category at time 1)
# ##   - later_category (IUCN category at time 2)
# 
# ## Example for Korean data:
raw2.1<- raw2 %>% select(no., 분류군, 국명,평가년도...5, 평가범주...6, 평가년도...8, 평가범주...9, 변동현황,멸종위기종) %>%
  mutate( taxa=분류군, speciesName=국명,
          early_year=평가년도...5, 
          early_category=평가범주...6,
          later_year=평가년도...8, 
          later_category=평가범주...9,
          updateCategory=변동현황
  ) %>%
  mutate(early_category = extract_category(early_category),
         later_category = extract_category(later_category),
         early_year=as.numeric(early_year), later_year= as.numeric(later_year)
  )  %>% 
  select(no., taxa, speciesName, early_year,early_category,later_year, later_category, updateCategory) %>%
  filter(early_year > 2000 & later_year > 2000)
# 
# ## Example for English data (categories already in standard format):
# # df <- data.frame(
# #   speciesName = c("Species A", "Species B", "Species C", ...),
# #   early_category = c("LC", "VU", "EN", ...),
# #   later_category = c("NT", "VU", "CR", ...)
# # )
# 
# 
## View sample size
categoryGivenSamples<-raw2.1 %>% select(taxa)  %>% table()
categoryGivenSamplesNotDD<-raw2.1 %>% filter(updateCategory != "신규" & !is.na(updateCategory) & 
                                               !is.na(early_category) & !is.na(later_category) &
                                               early_category !="DD" & later_category !="DD" ) %>%
  select(taxa)  %>% table()


## You do the filtering to narrow down specific taxa (example: plants & insects)
df_taxa <- subset(raw2.1, taxa %in% c("관속식물"))  # change taxa list as needed
df_taxa
df_taxa2<-df_taxa %>% filter(updateCategory != "신규" & !is.na(updateCategory) & 
                               !is.na(early_category) & !is.na(later_category) &
                               early_category !="DD" & later_category !="DD" )

# #### 2. Analyses ######
# ## Build t1/t2 vectors from your filtered df
# Your code should have speciesName,early_category,later_category in a dataframe format. 


v  <- prepare_status_vectors(df_taxa2)
t1 <- v$t1
t2 <- v$t2
# 
# ## Compute the full (truth) ΔRLI for the interval
truth <- delta_rli(t1, t2, mode = "overlap")
cat("Full ΔRLI:", truth$delta, "\n")
cat("N overlap:", truth$n_overlap, "\n")
# 
# ## Method 1: Bootstrap across sample sizes and plot (grid-based)
# ## Good for visualization, but resolution limited by step_n
bt <- bootstrap_srli(t1, t2, R = 50000, mode = "overlap", n_min = 0, step_n = 10, seed = 1023)
plot_srli_single_panel(bt, wrong_thresh = 0.05,
                       xlim=c(0,100),
                       ylim_left=c(0,20))


plot_srli_two_panel(bt, wrong_thresh = 0.05)
# 
# ## Grid-based recommendation (limited by step_n resolution)
n_rec_grid <- recommend_n_sign(bt, wrong_thresh = 0.05)
cat("Recommended n (grid-based):", n_rec_grid, "\n")
# 

# ## Method 2: Find exact threshold using binary search
# ## More accurate but takes longer; no visualization
result <- find_exact_n_threshold(t1, t2, wrong_thresh = 0.05, R = 50000, seed = 1023)
cat("Exact minimum n:", result$n_exact, "\n")
# cat("Error rate at n:", result$final_rate, "\n")


