In test1.opanal the analytic solution is a linear function.
With linear elements, the numerical solution is exact up to the accuracy
implied by the convergence tolerance of the iterative solver.

In test3.opanal the analytic solution is some sinusoidal function.
With only five linear elements, there is obviously a noticeable discretization
error. However, if beyond the 5^3 element mesh, we also compute the solution
on a 10^3 and 20^3 mesh, we see the expected convergence rates for velocity and pressure.
