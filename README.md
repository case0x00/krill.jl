# KRILL

(WIP)

this program simulates solutions to the circular restricted three body problem (CR3BP) in the synodic reference frame.

for more information about the theory, see [theory](##theory).


## zero velocity curves

using `plot_j_variation()`, plots variation of the jacobi integral for various ranges to demonstrate the effect it has on the zero velocity curves (hill region).

<p align="center"><img src="https://raw.githubusercontent.com/electric-coral/krill/master/src/CR3BP_jacobi_integral_variation.png"/></p>


## program functionality

* variation of the jacobi integral and animation of the change in zero velocity surface in the synodic frame of a 3 body system with lagrange points shown
* animation of a third body in the sidereal frame of a 3 body system for a given jacobi energy with lagrange points shown

## dependencies

* matplotlib
* tqdm

## usage

<pre>
$ ./run.sh
</pre>


## theory

the CR3BP is a variation of the three body problem (the dynamics of three gravitationally attracted bodies), whereby the orbits of the two larger bodies are circular and the third body has a negligable mass compared to the other two bodies. as such, the third body does not perturb the orbits of theother two bodies. hence circular (circular orbits) and restricted (third mass is restricted to be negligable).

the synodic reference frame is one whereby the two primary bodies lie stationary on the x-axis. in this example, the two body system is the earth-moon system which means the synodic frame rotates with the earth around the sun. in this frame, the third body's angular momentum and energy are not conserved, however there is a singular quantity that is conserved: the jacobi integral. variations of this value result in varying sizes of the hill region. a nice explanation is here at [space.stackexchange.com/q/21557](https://space.stackexchange.com/questions/21557/this-orbit-looks-wrong-near-a-lagrange-point-is-it/21570#21570).

the motion of the third body is restricted to regions whereby v^2 >= 0. these regions are bounded by surfaces of zero velocity. the value of the jacobi integral, J, determines the geometry of the ZVC (zero velocity curve). check out [this dissertation p40](https://engineering.purdue.edu/people/kathleen.howell.1/Publications/Dissertations/2011_CraigDavis.pdf) for a solid visualization of the variation of J with relation to J values at various lagrange points.



## resources

some good resources that I used

* [Craig Davis, D.E. 2011, MULTI-BODY TRAJECTORY DESIGN STRATEGIES BASED ON PERIAPSIS POINCARÉ MAPS -- PhD dissertation @ purdue](https://engineering.purdue.edu/people/kathleen.howell.1/Publications/Dissertations/2011_CraigDavis.pdf)
* [hyperphysics page on lagrange points](http://hyperphysics.phy-astr.gsu.edu/hbase/Mechanics/lagpt.html)
* [Westra, D. 2017 Lagrangian Points](https://www.mat.univie.ac.at/~westra/lagrangepoints.pdf)
* [nabla zero labs CR3BP sim](https://www.nablazerolabs.com/tbp/)
* [Chenciner, A. scholarpedia: three body problem](http://www.scholarpedia.org/article/Three_body_problem)
* [Meiss, J. scholarpedia: hamiltonian systems](http://www.scholarpedia.org/article/Hamiltonian_systems)
* [Gomez, D. and Barrabes, E. scholarpedia: space manifold dynamics](http://www.scholarpedia.org/article/Space_Manifold_dynamics)

