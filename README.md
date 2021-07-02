# N-body-problem

Solution of the N-Body Problem with 12-order Runge-Kutta in C using parallel processors with ALL-Pairs algorithm.

## Physics

Given $\mathrm{N}$ bodies in an inertial reference frame in $\mathbb{R}^{3}$ with an initial position $\overrightarrow{x_{i}}$ and velocity $\overrightarrow{v_{i}}$ for $1 \leq i \leq N$, the force vector $\overrightarrow{f_{i j}}$ on body $i$ caused by its gravitational attraction to body $j$ is given by the following equation:
$$
\overrightarrow{f_{i j}}=G \frac{m_{i} m_{j}}{\left|\vec{r}_{i j}\right|^{2}} \cdot \frac{\vec{r}_{i j}}{\left|\vec{r}_{i j}\right|}
$$
Where $m_{i}$ and $m_{j}$ are the masses of the bodies $i$ and $j$ respectively; and $\vec{r}_{i j}=\overrightarrow{x_{j}}-\vec{x}_{i}$ is the distance between the position vectors.
Summing over all masses yields the $N$ -body equations of motion:
$$
\vec{F}_{i}=m_{i} \frac{d^{2} \vec{x}_{i}}{d t^{2}}=\sum_{j=1 \atop j \neq i}^{n} \overrightarrow{f_{i j}}=\sum_{j=1 \atop j \neq i}^{n} \frac{G m_{i} m_{j}\left(\vec{r}_{i j}\right)}{\left\|\vec{r}_{i j}\right\|^{3}}=-\frac{\partial U}{\partial \vec{x}_{i}}
$$
where $U$ is the self-potential energy
$$
U=-\sum_{1 \leq i<j \leq n} \frac{G m_{i} m_{j}}{\left\|\vec{r}_{i j}\right\|}
$$


## Algorithm
The computation has two steps:
1. Compute the forces on each element
2. Move each element a bit based on this force, and then repeat

Since each object computes the forces on it from each other object ($N$ objects). Thus, the work complexity of this problem is a $N^2$ algorithm. Each object must compute its interaction with each other object so each of N objects has to compute $N−1$ forces.

<img src="img/N_body_problem_solar_system.jpg" width="750">