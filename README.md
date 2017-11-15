# SelectionExample

Worked examples of treating selection effects via sampling over
non-detections in Stan.  Contents:

 * Three different ways of hierarchically modelling a flux-limited
   survey of a population of objects with a log-normal luminosity
   function, with known photometric uncertainty in each observation.
   In this case the selection function can be computed analytically;
   methods that using sampling to account for the selection function
   are derived and validated by comparison to the analytic model.  See
   `Gaussian.ipynb`.

 * A more realistic worked example of inferring the parameters of a
   Schechter luminosity function from a hypothetical flux-limited
   galaxy survey in the nearby universe, again incorporating
   photometic uncertainties.  In this case, the selection function is
   not analytically tractible, so must be treated numerically.  See
   `Schechter.ipynb`.
