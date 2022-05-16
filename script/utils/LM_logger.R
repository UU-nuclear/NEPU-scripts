# This is a generator function that creates a logger function
# for the Levenberg-Marquardt optimization. This logger function
# can be passed as an argument to the LMalgo function.

createLoggerLM <- function(talys, saveDir)  {

    function(buf) {

      # this should be cached
      propDt <- talys$fun(buf$pprop, ret.dt = TRUE)
      setkey(propDt, IDX)
     
      # for this step it is important that 'talys'
      # caches calculation results to avoid 
      # unnecessary recalculations
      fref <- talys$fun(buf$pref, applySexp = FALSE)
      Jref <- talys$jac(buf$pref, applySexp = FALSE)
      
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



