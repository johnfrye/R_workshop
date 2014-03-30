GenLogistic <- function(time, A, K, b, v, Q, M) {
  # Compute the values of a generalized logistic function for a given set of 
  # times
  #
  # Args:
  #   time: Vector containing time values at which to compute the logistic 
  #     function.
  #   A: Lower asymptote of the logistic function.
  #   K: Upper asymptote of the logistic function.
  #   b: Growth rate of the logistic function. 
  #   v: Affects near which asymptote maximum growth occurs. v > 0.
  #   Q: Depends on the value at time 0.
  #   M: The time of maximum growth if Q=ν. 
  #
  # Returns:
  #   A vector of the same length as "time". 
  
  A + (K - A) / (1 + Q * exp(-b * (time - M))) ^ (1 / v)
}