# Historical biogeography analysis of Sapotaceae using BioGeoBEARS and nuclear data.
# For running with the plastome data, just substitute phylo and area files
#
# This script fits three models to a phylogeny
# and a set of geographic ranges, compares them statistically, and produces
# ancestral-range reconstructions:
#   1. DEC          
#   2. DIVALIKE     
#   3. BAYAREALIKE
#
# For each model the script:
#   - runs the maximum-likelihood optimization
#   - saves the raw results object (.Rdata)
#   - plots ancestral state reconstructions to PDF
#
# It then compares the three models using likelihood-ratio tests and AIC,
# and writes a summary table of model fit statistics to CSV/TXT.
#
# Usage:
#   Rscript biogeography.R <phylo.treefile> <areas.txt>
#
#
# Winter 2025
# RC Pizzardo
# =============================================================================

library(ape)          # phylogenetic tree handling
library(cladoRcpp)     # fast cladogenetic likelihood calculations (used internally by BioGeoBEARS)
library(BioGeoBEARS)   # core biogeographic model-fitting and plotting functions

# -----------------------------------------------------------------------------
# 1. Parse command-line arguments
# -----------------------------------------------------------------------------
argv <- commandArgs(trailingOnly = TRUE)
trfn   <- argv[1]  # path to the tree file
geogfn <- argv[2]  # path to the geographic range file

# -----------------------------------------------------------------------------
# 2. Read in the tree and geographic range data
# -----------------------------------------------------------------------------
moref(trfn)                    # print the tree file contents to the console
tr <- read.tree(trfn)          # load the tree as an ape "phylo" object
moref(geogfn)                  # print the geography file contents to the console
tipranges <- getranges_from_LagrangePHYLIP(lgdata_fn = geogfn)
max(rowSums(dfnums_to_numeric(tipranges@df)))
# Maximum number of areas any ancestral range is allowed to occupy.
max_range_size <- 5
# Report how many possible geographic states exist given 9 areas and a max
# range size of 5 (informational only - does not affect the analysis below)
numstates_from_numareas(numareas = 9, maxareas = 5, include_null_range = TRUE)

# =============================================================================
# 3. MODEL 1: DEC
# =============================================================================

# --- Build and configure the BioGeoBEARS run object -------------------------
BioGeoBEARS_run_object <- define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn              <- trfn
BioGeoBEARS_run_object$geogfn            <- geogfn
BioGeoBEARS_run_object$max_range_size    <- max_range_size
BioGeoBEARS_run_object$min_branchlength  <- 0.000001  # treat near-zero branches as this length
BioGeoBEARS_run_object$include_null_range <- TRUE     # allow the "no area occupied" state
BioGeoBEARS_run_object$on_NaN_error      <- -1e50      # heavily penalize NaN log-likelihoods during optimization
BioGeoBEARS_run_object$speedup           <- TRUE       # use faster (approximate) likelihood shortcuts
BioGeoBEARS_run_object$use_optimx        <- "GenSA"    # optimizer: Generalized Simulated Annealing
BioGeoBEARS_run_object$num_cores_to_use  <- 6          # parallel cores; adjust to match your machine
BioGeoBEARS_run_object$force_sparse      <- FALSE

# Load the tree and geography files into the run object
BioGeoBEARS_run_object <- readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Request extra outputs needed for downstream stats and ancestral state plots
BioGeoBEARS_run_object$return_condlikes_table               <- TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table <- TRUE
BioGeoBEARS_run_object$calc_ancprobs                         <- TRUE

# Inspect the run object and its default model parameters (printed to console)
BioGeoBEARS_run_object
BioGeoBEARS_run_object$BioGeoBEARS_model_object
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table

# Validate that inputs (tree, geography, parameters) are internally consistent
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

# --- Run the DEC analysis ----------------------------------------------------
resfn <- "sapotaceae_DEC_max5.Rdata"   # output file for the raw results object

res <- bears_optim_run(BioGeoBEARS_run_object)  # maximum-likelihood optimization
save(res, file = resfn)                          # persist raw results to disk
resDEC <- res                                    # keep a named copy for later comparison

# --- Plot the DEC ancestral state reconstruction -----------------------------
pdffn <- "sapotaceae_reconstruction_DEC_max5.pdf"
pdf(pdffn, width = 18, height = 45)

analysis_titletxt <- "BioGeoBEARS DEC on Sapotaceae"
results_object <- resDEC
scriptdir <- np(system.file("extdata/a_scripts", package = "BioGeoBEARS"))

# Reconstruction as text labels at each node
plot_BioGeoBEARS_results(results_object, analysis_titletxt, plotwhat = "text",
                          label.offset = 0.45, tipcex = 0.7, statecex = 0.7,
                          splitcex = 0.6, titlecex = 0.8, plotsplits = TRUE,
                          cornercoords_loc = scriptdir, include_null_range = TRUE,
                          tr = tr, tipranges = tipranges)

# Reconstruction as pie charts (state probabilities) at each node
plot_BioGeoBEARS_results(results_object, analysis_titletxt, plotwhat = "pie",
                          label.offset = 0.45, tipcex = 0.7, statecex = 0.7,
                          splitcex = 0.6, titlecex = 0.8, plotsplits = TRUE,
                          cornercoords_loc = scriptdir, include_null_range = TRUE,
                          tr = tr, tipranges = tipranges)

dev.off()

# =============================================================================
# 4. MODEL 2: DIVALIKE
# =============================================================================

# --- Build and configure the BioGeoBEARS run object -------------------------
BioGeoBEARS_run_object <- define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn              <- trfn
BioGeoBEARS_run_object$geogfn            <- geogfn
BioGeoBEARS_run_object$max_range_size    <- max_range_size
BioGeoBEARS_run_object$min_branchlength  <- 0.000001
BioGeoBEARS_run_object$include_null_range <- TRUE
BioGeoBEARS_run_object$on_NaN_error      <- -1e50
BioGeoBEARS_run_object$speedup           <- TRUE
BioGeoBEARS_run_object$use_optimx        <- "GenSA"
BioGeoBEARS_run_object$num_cores_to_use  <- 6          # You may have to change this depending on your computer
BioGeoBEARS_run_object$force_sparse      <- FALSE
# NOTE: these two lines override the settings above, switching to the
# built-in optimx optimizer and single-core execution.
BioGeoBEARS_run_object$use_optimx        <- TRUE
BioGeoBEARS_run_object$num_cores_to_use  <- 1
BioGeoBEARS_run_object$force_sparse      <- FALSE

BioGeoBEARS_run_object <- readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table               <- TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table <- TRUE
BioGeoBEARS_run_object$calc_ancprobs                         <- TRUE

# --- Constrain the model parameters to match the DIVALIKE model -------------
# DIVALIKE disallows "sympatric subset" speciation (s = 0)...
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s", "type"] <- "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s", "init"] <- 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s", "est"]  <- 0.0

# ...and splits cladogenetic weight equally between sympatry (y) and vicariance (v)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv", "type"] <- "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys",  "type"] <- "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y",   "type"] <- "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v",   "type"] <- "ysv*1/2"

# Allow classic, widespread vicariance; all cladogenetic events equiprobable
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v", "type"] <- "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v", "init"] <- 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v", "est"]  <- 0.5

# Enforce that fixed parameter values fall within their allowed min/max bounds
BioGeoBEARS_run_object <- fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object = BioGeoBEARS_run_object)
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

# --- Run the DIVALIKE analysis -----------------------------------------------
resfn <- "sapotaceae_DIVA_max5.Rdata"

res <- bears_optim_run(BioGeoBEARS_run_object)
save(res, file = resfn)
resDIVALIKE <- res

# --- Plot the DIVALIKE ancestral state reconstruction ------------------------
pdffn <- "sapotaceae_DIVA_max5.pdf"
pdf(pdffn, width = 18, height = 45)

analysis_titletxt <- "BioGeoBEARS DIVALIKE on Sapotaceae"
results_object <- resDIVALIKE
scriptdir <- np(system.file("extdata/a_scripts", package = "BioGeoBEARS"))

plot_BioGeoBEARS_results(results_object, analysis_titletxt, plotwhat = "text",
                          label.offset = 0.45, tipcex = 0.7, statecex = 0.7,
                          splitcex = 0.6, titlecex = 0.8, plotsplits = TRUE,
                          cornercoords_loc = scriptdir, include_null_range = TRUE,
                          tr = tr, tipranges = tipranges)

plot_BioGeoBEARS_results(results_object, analysis_titletxt, plotwhat = "pie",
                          label.offset = 0.45, tipcex = 0.7, statecex = 0.7,
                          splitcex = 0.6, titlecex = 0.8, plotsplits = TRUE,
                          cornercoords_loc = scriptdir, include_null_range = TRUE,
                          tr = tr, tipranges = tipranges)

dev.off()

# =============================================================================
# 5. MODEL 3: BAYAREALIKE
# =============================================================================

# --- Build and configure the BioGeoBEARS run object -------------------------
BioGeoBEARS_run_object <- define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn              <- trfn
BioGeoBEARS_run_object$geogfn            <- geogfn
BioGeoBEARS_run_object$max_range_size    <- max_range_size
BioGeoBEARS_run_object$min_branchlength  <- 0.000001
BioGeoBEARS_run_object$include_null_range <- TRUE
BioGeoBEARS_run_object$on_NaN_error      <- -1e50
BioGeoBEARS_run_object$speedup           <- TRUE
BioGeoBEARS_run_object$use_optimx        <- "GenSA"
BioGeoBEARS_run_object$num_cores_to_use  <- 6          # You may have to change this depending on your computer
BioGeoBEARS_run_object$force_sparse      <- FALSE
# NOTE: as above, these overrides switch to optimx / single-core execution.
BioGeoBEARS_run_object$use_optimx        <- TRUE
BioGeoBEARS_run_object$num_cores_to_use  <- 1
BioGeoBEARS_run_object$force_sparse      <- FALSE

BioGeoBEARS_run_object <- readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table               <- TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table <- TRUE
BioGeoBEARS_run_object$calc_ancprobs                         <- TRUE

# --- Constrain the model parameters to match the BAYAREALIKE model ----------
# No "sympatric subset" speciation (s = 0)...
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s", "type"] <- "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s", "init"] <- 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s", "est"]  <- 0.0

# ...and no vicariance (v = 0): BayArea assumes ranges are simply copied at speciation
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v", "type"] <- "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v", "init"] <- 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v", "est"]  <- 0.0

# Adjust linkage between the remaining cladogenetic parameters accordingly
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv", "type"] <- "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys",  "type"] <- "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y",   "type"] <- "1-j"

# Only sympatric/range-copying (y) events allowed, with near-exact copying
# (both descendants essentially always inherit the full ancestral range)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y", "type"] <- "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y", "init"] <- 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y", "est"]  <- 0.9999

# Check the inputs; fixing any initial ("init") values outside their min/max bounds
BioGeoBEARS_run_object <- fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object = BioGeoBEARS_run_object)
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

# --- Run the BAYAREALIKE analysis --------------------------------------------
resfn <- "sapotaceae_BAYAREALIKE_max5.Rdata"

res <- bears_optim_run(BioGeoBEARS_run_object)
save(res, file = resfn)
resBAYAREALIKE_null <- res

# --- Plot the BAYAREALIKE ancestral state reconstruction ---------------------
pdffn <- "sapotaceae_BAYAREALIKE_max5.pdf"
pdf(pdffn, width = 18, height = 45)

analysis_titletxt <- "BioGeoBEARS BAYAREALIKE on Sapotaceae"
results_object <- resBAYAREALIKE_null
scriptdir <- np(system.file("extdata/a_scripts", package = "BioGeoBEARS"))

plot_BioGeoBEARS_results(results_object, analysis_titletxt, plotwhat = "text",
                          label.offset = 0.45, tipcex = 0.7, statecex = 0.7,
                          splitcex = 0.6, titlecex = 0.8, plotsplits = TRUE,
                          cornercoords_loc = scriptdir, include_null_range = TRUE,
                          tr = tr, tipranges = tipranges)

plot_BioGeoBEARS_results(results_object, analysis_titletxt, plotwhat = "pie",
                          label.offset = 0.45, tipcex = 0.7, statecex = 0.7,
                          splitcex = 0.6, titlecex = 0.8, plotsplits = TRUE,
                          cornercoords_loc = scriptdir, include_null_range = TRUE,
                          tr = tr, tipranges = tipranges)

dev.off()

# =============================================================================
# 6. Model comparison statistics
# =============================================================================
# Pairwise likelihood-ratio-style comparisons between the three models.
# Each pair differs by exactly one free parameter ("j", the founder-event
# speciation parameter), hence numparams = 2 for each model pair below.

# --- Extract log-likelihoods --------------------------------------------------
LnL_1 <- get_LnL_from_BioGeoBEARS_results_object(resDIVALIKE)
LnL_2 <- get_LnL_from_BioGeoBEARS_results_object(resDEC)
LnL_3 <- get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKE_null)

# Number of free parameters estimated in each model (here: dispersal "d",
# extinction "e" - or their DIVALIKE/BAYAREALIKE equivalents)
numparams1 <- 2  # DIVALIKE
numparams2 <- 2  # DEC
numparams3 <- 2  # BAYAREALIKE

# Pairwise AIC/LRT-style comparisons
stats  <- AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)  # DIVALIKE vs DEC
stats2 <- AICstats_2models(LnL_3, LnL_2, numparams3, numparams2)  # BAYAREALIKE vs DEC
stats3 <- AICstats_2models(LnL_3, LnL_1, numparams3, numparams1)  # BAYAREALIKE vs DIVALIKE

# --- Extract per-model parameter estimates (d, e, j, LnL, etc.) -------------
res2 <- extract_params_from_BioGeoBEARS_results_object(results_object = resDEC,
          returnwhat = "table", addl_params = c("j"), paramsstr_digits = 4)
res1 <- extract_params_from_BioGeoBEARS_results_object(results_object = resDIVALIKE,
          returnwhat = "table", addl_params = c("j"), paramsstr_digits = 4)
res3 <- extract_params_from_BioGeoBEARS_results_object(results_object = resBAYAREALIKE_null,
          returnwhat = "table", addl_params = c("j"), paramsstr_digits = 4)

# Quick console preview of the combined parameter tables
rbind(res2, res1)
rbind(res2, res3)

# Format the LRT-style test statistics for reporting
tmp_tests  <- conditional_format_table(stats)
tmp_tests2 <- conditional_format_table(stats2)

# --- Assemble summary tables --------------------------------------------------
restable    <- NULL
teststable  <- NULL
restable2   <- NULL
teststable2 <- NULL

# Table 1: DEC vs DIVALIKE
restable   <- rbind(restable, res2, res1)
teststable <- rbind(teststable, tmp_tests)
teststable$alt  <- c("DIVALIKE")
teststable$null <- c("DEC")
row.names(restable) <- c("DEC", "DIVALIKE")
restable <- put_jcol_after_ecol(restable)  # reorder columns so "j" appears after "e"

# Table 2: DEC vs BAYAREALIKE
restable2   <- rbind(restable2, res2, res3)
teststable2 <- rbind(teststable2, tmp_tests2)
teststable2$alt  <- c("BAYAREALIKE")
teststable2$null <- c("DEC")
row.names(restable2) <- c("DEC", "BAYAREALIKE")
restable2 <- put_jcol_after_ecol(restable2)

# --- Add AIC and Akaike weights ----------------------------------------------
AICtable <- calc_AIC_column(LnL_vals = restable$LnL, nparam_vals = restable$numparams)
restable <- cbind(restable, AICtable)
restable_AIC_rellike <- AkaikeWeights_on_summary_table(restable = restable, colname_to_use = "AIC")
restable_AIC_rellike <- put_jcol_after_ecol(restable_AIC_rellike)

AICtable2 <- calc_AIC_column(LnL_vals = restable2$LnL, nparam_vals = restable2$numparams)
restable2 <- cbind(restable2, AICtable2)
restable_AIC_rellike2 <- AkaikeWeights_on_summary_table(restable = restable2, colname_to_use = "AIC")
restable_AIC_rellike2 <- put_jcol_after_ecol(restable_AIC_rellike2)

# --- Combine and save the final model comparison table ----------------------
# Combines DEC, DIVALIKE, and BAYAREALIKE into one table (DEC's row from the
# first table is kept; only BAYAREALIKE's row is appended from the second).
stats_table <- rbind(restable_AIC_rellike, restable_AIC_rellike2[2, ])
stats_table

write.csv(as.data.frame(stats_table), "restable_AIC_rellike.csv")
write.table(conditional_format_table(stats_table), file = "restable_AIC_rellike.txt",
            sep = "\t", quote = FALSE)

### THE END ###
