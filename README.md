# SelectionExample

Worked examples of treating selection effects, including via sampling over
non-detections in Stan.  See [Mandel, Farr & Gair
(2018)](https://ui.adsabs.harvard.edu/#abs/2018arXiv180902063M/abstract).
Contents:

 * A simple worked example of a flux-limited survey of standard candles in
   Euclidean space illustrating that the commonly-used modification of the
   population distribution to incorporate the selection function generates
   biased inferences over population parameters.  The correct solution is
   briefly discussed and demonstrated to be correct.

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
