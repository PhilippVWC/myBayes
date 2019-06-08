#  Logit-Normal-Verteilung
#  Wird parametrisiert durch die Parameter mu und sd der zugrundeliegenden Normalverteilung
#  Fuer die Momente der Logit-Normal-Verteilung gibt es allerdings keine analytische Loesung.
#  Gesucht waere ein Funktion, die aus vorgegebenem Mittelwert und SD,
#  die entsprechenden mu, sd der Normalverteilung berechnet.

#' @export
rlogitnormal = function(...){
  x = exp(rnorm(...))
  x / (1+x)
}

#samp = rlogitnormal(1000, 0, 1)
#summary(samp);1/(1+exp(0))
#samp = rlogitnormal(1000, 1, 1)
#summary(samp);1/(1+exp(-1))
#samp = rlogitnormal(1000, 1, 10)
#summary(samp);1/(1+exp(-1))
