LogBack <- function(x, w, alph, pfcn) {
  return(pfcn(x) * x^alph * (w - x)^alph / w^(2 * alph + 1))
}
