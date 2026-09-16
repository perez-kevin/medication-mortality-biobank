library(survival)
library(progressr)
library(survival)
library(broom)
library(dplyr)
library(purrr)

# Example data.frame
df <- data.frame(
  ID           = c("p1","p2","p3","p3","p4","p5","p6"),
  Gender       = c("F","M","F","F","M","F","M"),
  YearEnrolled = c(2004,2006,2005,2005,2004,2007,2006),
  AgeEnroll    = c(61,55,72,72,48,66,59),
  Time         = c(8.2,5.0,3.4,3.4,10.0,6.7,2.1),   # follow-up, years
  Status       = c(1,0,1,1,0,1,1),                   # 1 = died
  Normalized   = c("aspirin","aspirin","aspirin","statin","statin","statin","aspirin"),
  TimeToDrug   = c(2.0,1.5,0.8,2.9,4.0,3.3,0.5)      # years from enrollment to first script
)

run_subset_fast <- function(subset_df, drug_list, group_name, progressor_fn, formula) {
  
  patient_baseline <- subset_df %>%
    distinct(ID, Gender, YearEnrolled, AgeEnroll, Time, Status, .keep_all = T)
  
  if(length(unique(patient_baseline$Gender)) > 1) {
    f <- as.formula(formula)
  } else {
    f <- as.formula("Surv(tstart, tstop, status) ~ using + AgeEnroll")
  }
  
  map_dfr(drug_list, function(drug) {
    progressor_fn()
    
    this_drug_events <- subset_df %>%
      filter(Normalized == drug) %>%
      select(ID, TimeToDrug)
    
    if (nrow(this_drug_events) == 0)
      return(NULL)
    
    model_data <- tmerge(patient_baseline, patient_baseline, id = ID, 
                         status = event(Time, Status))
    
    model_data <- tmerge(model_data, this_drug_events, id = ID, 
                         using = tdc(TimeToDrug))
    
    MIN_DEATHS = 0 # for example df
    # MIN_DEATHS = 15 # toggle with real data
    
    if (sum(model_data$using & model_data$Status) < MIN_DEATHS)
      return(NULL)
    
    model <- coxph(f, data = model_data)
    res <- tidy(model, conf.int = TRUE) %>% filter(term == "using")
    
    tibble(
      drug    = drug,
      HR      = exp(res$estimate),
      conf.low  = exp(res$conf.low),
      conf.high = exp(res$conf.high),
      p.value = res$p.value,
      Group   = group_name
    )
  })
}

handlers(handler_progress(
  format = "[:bar] :percent :eta | :message",
  clear = FALSE, 
  width = 60
))


scenarios <- list(
  "Total"  = df
  # "Diabetes" = df %>% filter(Diabetes),
  # "Cancer"   = df %>% filter(Cancer),
  # "CVD"   = df %>% filter(CVD),
  # "Smoker"   = df %>% filter(Smoker),
  # "Healthy"   = df %>% filter(!(CVD | Cancer | Diabetes)),
  # "Male"   = df %>% filter(Gender == 'Male'),
  # "Female"   = df %>% filter(Gender == 'Female')
)

# 1. Define your list of formula strings
formulas <- list(
  "Model1" = "Surv(tstart, tstop, status) ~ using + AgeEnroll + Gender + YearEnrolled"
  # "Model2" = "Surv(tstart, tstop, status)  ~ using + AgeEnroll + Gender + YearEnrolled + CCI",
  # "Model3" = "Surv(tstart, tstop, status)  ~ using + AgeEnroll + Gender + YearEnrolled + CCI + Smoker + BMI + SES"
)

drug_list = c('statin', 'aspirin')

# 2. Run the nested iteration
all_results <- imap_dfr(scenarios, function(data, scen_name) {
  
  # Iterate over each formula for the current scenario
  imap_dfr(formulas, function(formula_str, form_name) {
    
    # Optional: create a descriptive label for the progress bar
    current_label <- paste0(scen_name, " (", form_name, ")")
    
    # Wrap the specific run in progress handling
    results <- with_progress({
      p <- progressor(steps = length(drug_list), message = sprintf("Running %s", current_label))
      
      # Pass the formula_str to your function
      run_subset_fast(data, drug_list, scen_name, p, formula = formula_str)
    })
    
    # print(results)
    
    # Add metadata columns and adjust p-values
    results %>%
      mutate(
        scenario = scen_name,
        model_type = form_name,
        formula_used = formula_str,
        # Adjust p-values based on the number of tests in THIS specific subset
        p.adj = pmin(p.value * n(), 1)
      )
  })
}) %>%
  # Relocate new columns to the front for better readability
  relocate(scenario, model_type, p.adj, .after = p.value)

all_results$Group = all_results$model_type
# all_results$Group = all_results$scenario
all_results$model_type = NULL; all_results$scenario = NULL; all_results$formula_used = NULL
all_results
