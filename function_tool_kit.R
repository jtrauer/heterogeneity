
# stand-alone functions for use by other scripts

# simple function to convert the LHS sampling over the (0, 1) interval to a sample
# from the range that is specific to the particular parameter
adjust_lhs_to_range <- function(input_value, parameter_name, parameter_ranges) {
    max <- parameter_ranges[[parameter_name]]$max
    min <- parameter_ranges[[parameter_name]]$min
    output_value <- input_value * (max - min) + min
}

# find the flow rate of delta based on a case detection rate (proportion), other outflows
# from the compartment and a final adjustment for the treatment success rate (proportion)
find_delta_from_cdr <- function(cdr, outflows, treatment_success) {
  cdr * outflows / (1 - cdr) * treatment_success
}

# very simple function to get rate from proportion and sojourn time
# pass parameters to similar function to simplify input code
find_rate_from_proportion_params <- function(proportion, sojourn_time, parameters, modifier="") {
  if (modifier == "") {
    modifier <- 1
    }
  else {
    modifier <- parameters[[modifier]]
    }
  parameters[[proportion]] / parameters[[sojourn_time]] * modifier
}

