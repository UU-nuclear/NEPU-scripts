createLoggerLMalt <- function(saveDir) {

	function(buf) {
      cat("### iteration: ", buf$iteration, "\n",
          "mu:        ", buf$mu, "\n",
          "gain:      ", buf$gain, "\n",
          "Lref:      ", buf$Lref, "\n",
          "Lprop:     ", buf$Lprop, "\n",
          "Lprop_apx: ", buf$Lprop_approx, "\n",
          file = file.path(saveDir, "lastLM.log"), append = TRUE)
      cat("### iteration:", buf$iteration, "\n",
          "pars: ", buf$pprop, "\n",
          file = file.path(saveDir, "lastPars.log"), append = TRUE)
    }
    
}