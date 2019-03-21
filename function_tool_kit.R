
# stand-alone functions for use by other scripts

adjust_lhs_to_range <- function(input_value, parameter_name, parameter_ranges) {
    max <- parameter_ranges[[parameter_name]]$max
    min <- parameter_ranges[[parameter_name]]$min
    output_value <- input_value * (max - min) + min
}
