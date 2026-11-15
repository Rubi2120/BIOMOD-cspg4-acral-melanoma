# BIOMOD 2026 - Pandora's Box
# Exploratory overall survival analysis of CSPG4 expression
# in 44 patients with primary acral melanoma.
#
# This script reproduces the results presented at BIOMOD.
# Survival is measured from diagnosis, without delayed entry.
#
# Clinical dates have month-year resolution.
# Three-digit years are interpreted by appending a trailing
# zero, following the original analysis. This assumption
# remains to be confirmed.
#
# Individual-level input data are not publicly distributed.

library(survival)

input_root <- getOption(
  "biomod.survival.input_root", "."
)

output_dir <- getOption(
  "biomod.survival.output_dir", "tables"
)

# Load original expression and clinical data

expression <- read.csv(
  file.path(
    input_root,
    "results_v2/tables/nivelA_CSPG4_TPMcrudo_86muestras.csv"
  ),
  check.names = FALSE
)

expression <- expression[expression$en_44, ]

clinical <- read.csv(
  file.path(
    input_root,
    "data/raw/Supplementary_Table_1.csv"
  ),
  check.names = FALSE
)

# Match expression and clinical records

expression$Sample <- sub(
  "^PR", "PD", expression$sample_id
)

dat <- merge(
  expression[, c("Sample", "CSPG4_log2")],
  clinical[, c(
    "Sample",
    "Patient",
    "Date_of_diagnosis",
    "Date_of_last_note",
    "Date_of_death"
  )],
  by = "Sample"
)

cat(
  "Matched samples:", nrow(dat),
  "of 44 | Unique patients:",
  length(unique(dat$Patient)), "\n"
)

if (anyDuplicated(dat$Patient)) {
  dat <- dat[!duplicated(dat$Patient), ]
  cat("One sample per patient:", nrow(dat), "\n")
}

# Convert month-year dates to absolute months

parse_month <- function(x) {

  x <- trimws(as.character(x))
  valid <- !is.na(x) & x != ""
  result <- rep(NA_real_, length(x))

  parts <- strsplit(
    x[valid], ".", fixed = TRUE
  )

  month <- sapply(
    parts, function(z) as.numeric(z[1])
  )

  year <- sapply(
    parts,
    function(z) {
      as.numeric(substr(
        paste0(z[2], "000"), 1, 4
      ))
    }
  )

  stopifnot(
    all(month >= 1 & month <= 12),
    all(year >= 2000 & year <= 2026)
  )

  result[valid] <- year * 12 + month

  result
}

# Define overall survival time and event status

diagnosis <- parse_month(dat$Date_of_diagnosis)
last_note <- parse_month(dat$Date_of_last_note)
death <- parse_month(dat$Date_of_death)

dat$event <- as.integer(!is.na(death))

dat$time <- ifelse(
  !is.na(death), death, last_note
) - diagnosis

cat(
  "Positive survival time:",
  sum(dat$time > 0, na.rm = TRUE),
  "| Death events:",
  sum(dat$event[dat$time > 0], na.rm = TRUE),
  "\n"
)

dat <- dat[
  !is.na(dat$time) & dat$time > 0,
]

# Continuous-expression Cox model

cox_model <- coxph(
  Surv(time, event) ~ CSPG4_log2,
  data = dat
)

model_summary <- summary(cox_model)

# Median-expression groups and log-rank test

dat$group <- ifelse(
  dat$CSPG4_log2 > median(dat$CSPG4_log2),
  "high",
  "low"
)

logrank <- survdiff(
  Surv(time, event) ~ group,
  data = dat
)

# Export aggregate results only

results <- data.frame(
  n = nrow(dat),
  events = sum(dat$event),
  median_followup_months = median(dat$time),
  HR_per_log2_unit = model_summary$coefficients[1, 2],
  CI95_lower = model_summary$conf.int[1, 3],
  CI95_upper = model_summary$conf.int[1, 4],
  p_cox = model_summary$coefficients[1, 5],
  p_logrank_high_vs_low =
    1 - pchisq(logrank$chisq, 1),
  p_proportional_hazards = tryCatch(
    cox.zph(cox_model)$table[1, 3],
    error = function(e) NA
  )
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  results,
  file.path(
    output_dir,
    "acral_survival_biomod.csv"
  ),
  row.names = FALSE
)

print(results, digits = 5)

cat(
  "HR > 1 indicates an association between higher",
  "CSPG4 expression and a higher hazard of death.\n",
  "This exploratory model is unadjusted for age",
  "and clinical stage.\n"
)
