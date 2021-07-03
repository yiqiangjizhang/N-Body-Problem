# N-body-problem

Solution of the N-Body Problem with 12-order Runge-Kutta in C using parallel processors with ALL-Pairs algorithm.

## Physics

Given <!-- $\mathrm{N}$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cmathrm%7BN%7D"> bodies in an inertial reference frame in <!-- $\mathbb{R}^{3}$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cmathbb%7BR%7D%5E%7B3%7D"> with an initial position <!-- $\overrightarrow{x_{i}}$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Coverrightarrow%7Bx_%7Bi%7D%7D"> and velocity <!-- $\overrightarrow{v_{i}}$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Coverrightarrow%7Bv_%7Bi%7D%7D"> for <!-- $1 \leq i \leq N$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=1%20%5Cleq%20i%20%5Cleq%20N">, the force vector <!-- $\overrightarrow{f_{i j}}$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Coverrightarrow%7Bf_%7Bi%20j%7D%7D"> on body <!-- $i$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=i"> caused by its gravitational attraction to body <!-- $j$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=j"> is given by the following equation:
<!-- $$
\overrightarrow{f_{i j}}=G \frac{m_{i} m_{j}}{\left|\vec{r}_{i j}\right|^{2}} \cdot \frac{\vec{r}_{i j}}{\left|\vec{r}_{i j}\right|}
$$ --> 

<div align="center"><img style="background: white;" src="https://render.githubusercontent.com/render/math?math=%5Coverrightarrow%7Bf_%7Bi%20j%7D%7D%3DG%20%5Cfrac%7Bm_%7Bi%7D%20m_%7Bj%7D%7D%7B%5Cleft%7C%5Cvec%7Br%7D_%7Bi%20j%7D%5Cright%7C%5E%7B2%7D%7D%20%5Ccdot%20%5Cfrac%7B%5Cvec%7Br%7D_%7Bi%20j%7D%7D%7B%5Cleft%7C%5Cvec%7Br%7D_%7Bi%20j%7D%5Cright%7C%7D%0D"></div>

Where <!-- $m_{i}$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=m_%7Bi%7D"> and <!-- $m_{j}$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=m_%7Bj%7D"> are the masses of the bodies <!-- $i$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=i"> and <!-- $j$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=j"> respectively; and <!-- $\vec{r}_{i j}=\overrightarrow{x_{j}}-\vec{x}_{i}$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cvec%7Br%7D_%7Bi%20j%7D%3D%5Coverrightarrow%7Bx_%7Bj%7D%7D-%5Cvec%7Bx%7D_%7Bi%7D"> is the distance between the position vectors.
Summing over all masses yields the <!-- $N$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=N"> -body equations of motion:
<!-- $$
\vec{F}_{i}=m_{i} \frac{d^{2} \vec{x}_{i}}{d t^{2}}=\sum_{j=1 \atop j \neq i}^{n} \overrightarrow{f_{i j}}=\sum_{j=1 \atop j \neq i}^{n} \frac{G m_{i} m_{j}\left(\vec{r}_{i j}\right)}{\left\|\vec{r}_{i j}\right\|^{3}}=-\frac{\partial U}{\partial \vec{x}_{i}}
$$ --> 

<div align="center"><img style="background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cvec%7BF%7D_%7Bi%7D%3Dm_%7Bi%7D%20%5Cfrac%7Bd%5E%7B2%7D%20%5Cvec%7Bx%7D_%7Bi%7D%7D%7Bd%20t%5E%7B2%7D%7D%3D%5Csum_%7Bj%3D1%20%5Catop%20j%20%5Cneq%20i%7D%5E%7Bn%7D%20%5Coverrightarrow%7Bf_%7Bi%20j%7D%7D%3D%5Csum_%7Bj%3D1%20%5Catop%20j%20%5Cneq%20i%7D%5E%7Bn%7D%20%5Cfrac%7BG%20m_%7Bi%7D%20m_%7Bj%7D%5Cleft(%5Cvec%7Br%7D_%7Bi%20j%7D%5Cright)%7D%7B%5Cleft%5C%7C%5Cvec%7Br%7D_%7Bi%20j%7D%5Cright%5C%7C%5E%7B3%7D%7D%3D-%5Cfrac%7B%5Cpartial%20U%7D%7B%5Cpartial%20%5Cvec%7Bx%7D_%7Bi%7D%7D%0D"></div>

where <!-- $U$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=U"> is the self-potential energy
<!-- $$
U=-\sum_{1 \leq i<j \leq n} \frac{G m_{i} m_{j}}{\left\|\vec{r}_{i j}\right\|}
$$ --> 

<div align="center"><img style="background: white;" src="https://render.githubusercontent.com/render/math?math=U%3D-%5Csum_%7B1%20%5Cleq%20i%3Cj%20%5Cleq%20n%7D%20%5Cfrac%7BG%20m_%7Bi%7D%20m_%7Bj%7D%7D%7B%5Cleft%5C%7C%5Cvec%7Br%7D_%7Bi%20j%7D%5Cright%5C%7C%7D%0D"></div>


## Algorithm
The computation has two steps:
1. Compute the forces on each element
2. Move each element a bit based on this force, and then repeat

Since each object computes the forces on it from each other object (<!-- $N$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=N"> objects). Thus, the work complexity of this problem is a <!-- $N^2$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=N%5E2"> algorithm. Each object must compute its interaction with each other object so each of <!-- $N$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=N"> objects has to compute <!-- $N-1$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=N-1"> forces.

<img src="img/N_body_problem_solar_system.jpg" width="750">