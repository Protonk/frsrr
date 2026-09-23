# Independent conversions through R's binary I/O; no package kernel is used.
float32 <- function(x) readBin(writeBin(as.double(x), raw(), size = 4),
                               'double', n = length(x), size = 4)
float_from_bits <- function(bits) readBin(writeBin(as.integer(bits), raw(), size = 4),
                                         'double', n = length(bits), size = 4)
search_values <- function(x, magics, NRmax = 1L,
                          objective = 'max_relative_error', dependent = 'avg_relative_error') {
    .Call('_frsrr_search_optimal_constant', PACKAGE = 'frsrr', as.double(x),
          as.integer(magics), as.integer(NRmax), objective, dependent)
}
