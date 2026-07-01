# Prior and posterior predictive checks for the F2 food-group diversity model (E1).
#
# These were the two slow, redundant model runs in 172: it refit a `sample_prior =
# "only"` model for the prior check and then refit the main model again for the
# posterior check. Here both checks read fits already cached by 10:
#   - the prior predictive check uses the cached prior-only fit
#   - the posterior predictive check reuses the SAME main food-group fit that F2 is
#     drawn from (no extra refit), which is the whole point of caching it
# So nothing is sampled in this script; it only draws from stored fits.

source(here::here("reproduce", "human", "_human_helpers.R"))

prior_fit <- readRDS(cache_path("172_fit_fg_prior_only.rds"))
fg_fit    <- readRDS(cache_path("172_fit_fg_diversity.rds"))

# A bottom legend, stacked over two rows with small text, so the long
# "observed / simulated ..." labels are never clipped at the figure edge.
legend_two_rows <- function() {
  list(theme(legend.position = "bottom", legend.title = element_blank(),
             legend.text = element_text(size = 8),
             legend.key.size = unit(0.9, "lines")),
       guides(colour = guide_legend(nrow = 2, byrow = TRUE),
              fill   = guide_legend(nrow = 2, byrow = TRUE)))
}

# E1g: prior predictive check -------------------------------------------------
set.seed(8)
prior_pred <- pp_check(prior_fit, ndraws = 500, alpha = 0.1) +
  labs(y = "log(microbiome alpha diversity)", title = "Prior predictive check") +
  scale_color_discrete(labels = c("observed distribution",
                                  "simulated from prior predictive distribution")) +
  legend_two_rows()
save_panel(prior_pred, "E1g_prior_predictive_check.pdf", width = 170, height = 125)

# E1h: posterior predictive check (reuses the cached main fit) -----------------
# Published colour scheme: observed drawn in salmon/coral, posterior draws in
# teal, both with white box fills so only the outlines carry the colour (matched
# scale_color/scale_fill so the two legends merge into one).
post_pred <- pp_check(fg_fit, type = "boxplot", ndraws = 10, notch = FALSE) +
  labs(y = "log(microbiome alpha diversity)", title = "Posterior predictive check") +
  scale_color_manual(values = c("#F08A80", "#40B7AD"),
                     labels = c("observed distribution",
                                "simulated from posterior predictive distribution")) +
  scale_fill_manual(values = c("white", "white"),
                    labels = c("observed distribution",
                               "simulated from posterior predictive distribution")) +
  legend_two_rows()
save_panel(post_pred, "E1h_posterior_predictive_check.pdf", width = 170, height = 125)

message("Predictive-check diagnostics written to results/.")
