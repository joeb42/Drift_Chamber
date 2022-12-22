## The Drift Chamber

The drift chamber is an example of a wire chamber that have been used in high energy physics since the 1960s. Specifically it is used to detect ionising radiation and recontruct the paths of incident high energy radiation. Here the performance of a drift chamber is simulated in the context of the detection of cosmics rays. At altitudes ~ sea level the majority of cosmic rays are muons and so here it is assumed that all incident radiation are muon particles. 

### Basic diagram showing the set-up simulated (side view)

<div>
  <img width="400" alt="Screenshot 2021-05-15 at 16 26 02" src="https://user-images.githubusercontent.com/70596457/118369048-c8008380-b59a-11eb-93cf-5900b127d6aa.png">
  </div>

### Example of initial ionisation pattern

<div>
  
![Example_ion](https://user-images.githubusercontent.com/70596457/118369930-452cf800-b59d-11eb-8909-4f1cc50d8288.png)
</div>


## Monte Carlo Treatment of Cosmic Ray Muons

Incident cosmic ray muons were simulated in a monte carlo fashion, predominantly making use of the random function from the Numpy random library. Here we are interested in mapping the uniform probability distribution on the interval \[0, 1) provided by the random function to a general PDF P(y) on some specified interval between y<sub>0</sub> and y<sub>1</sub>. In order for this to work we require that the CDFs (cumulative density functions) are equal. This means that for some randomly generated x<sub>in</sub> from the random function we have the following equality <br>
![equation](https://latex.codecogs.com/gif.latex?%5Cint_%7Bx_0%7D%5E%7Bx_%7Bin%7D%7D%20u%28x%29%20dx%20%3D%20%5Cint_%7By_0%7D%5E%7By_%7Bout%7D%7D%20P%28y%29%20dy)
 <br> where u(x) is the uniform distribution provided by the random function and y<sub>out</sub> is the desired output. The LHS of this comes out as just x<sub>in</sub> which leads to the following equation for y<sub>out</sub> <br>
![equation](https://latex.codecogs.com/gif.latex?y_%7Bout%7D%20%3D%20F%5E%7B-1%7D%28x_%7Bin%7D%29)
<br>, where the function F is the RHS of the equality above. 

### Accept Reject Method

Of course this method is only valid where the function F is analytically invertable - where this is not the case we need to be a bit more creative. For example the zenith angle (theta) of the muon direction (modelled in spherical coordinates) is proportional to cos<sup>2</sup>(theta) which is not analytically invertable. Here we require what is known as an accept/reject method. Here we generate a pair of random values x1 and x2, the first of which is uniform on the domain of our pdf (so for the zenith angle this is \[pi/2, 3pi/2]) and the second uniform on the range of our pdf (\[0,1)). If x2 < P(x1) we 'accept' x1 and take it as our random value. Otherwise we continue the loop until an accepted pair is found. 

## Numerical simulation of ion drift diffusion 

After the Argon gas is ionised by the muon, the newly charged particles will undergo both drift (owing to the large electric field) and diffusion. Here the scenario is modelled numerically by discretising the chamber into a series of grid squares and modelling the changes in the charge distribution over time with a finite difference method. 
<br>
The change in the charge distribution q(x,y) over time is modelled by the following differential equation:
<br>
![equation](https://latex.codecogs.com/gif.latex?%5Cfrac%7B%5Cpartial%20q%7D%7B%5Cpartial%20t%7D%20%3D%20D%20%5Cnabla%5E2q%20-%20%5Cfrac%7B1%7D%7B%5Cmu%7D%20E%5Cnabla%20q)
<br> 
, where D is the diffusivity of the Argon gas that quantifies the amount of diffusion taking place. 
Using an implicit finite difference method one arrives at the following expression for the charge at grid coordinate (i,j) at time t:
<br>
![equation](https://latex.codecogs.com/gif.latex?q_%7Bi%2Cj%7D%5Et%20%3D%20%28-a-b%29q_%7Bi-1%2C%20j%7D%5E%7Bt&plus;1%7D%20-%20aq_%7Bi&plus;1%2Cj%7D%5E%7Bt&plus;1%7D%20-%20aq_%7Bi%2C%20j&plus;1%7D%5E%7Bt&plus;1%7D%20-%20aq_%7Bi%2C%20j-1%7D%5E%7Bt&plus;1%7D%20&plus;%20%281&plus;4a&plus;b%29q_%7Bi%2Cj%7D%5E%7Bt&plus;1%7D)
<br>
Here a = (D * timestep)/spacing and b = timestep * E / (2 * mu * spacing). 
This can be written compactly as the matrix equation
<br>
![equation](https://latex.codecogs.com/gif.latex?%5Cvec%7Bq%5Et%7D%20%3D%20M%5Cvec%7Bq%5E%7Bt&plus;1%7D%7D)
<br>
It is then a mere matter of solving a matrix equation at each time step. This can be optimised nicely by using SciPy's Sparse library, although memory considerations mean that care needs to be taken in not setting too small of a grid spacing. At present the primary limitation of this method is in the existence of a periodic boundary condition. That is, all charge that drift diffuses through a boundary reappears through the opposite boundary. Here the effect is damped by manually setting all charge to 0 within 2cm of the chamber but this is somewhat unsatifactory and needs further work.

## Running
Note this project depends on tkinter so make sure you are running a python distribution that is compatible (I built this project using miniconda). <br> 
If using miniconda ensure you have pip installed into your environment with ```conda install pip```. <br>
Run the command ```pip install -e .``` to install all dependencies. <br>
Run ```python main.py```: this should pop up a GUI window and you should be able to select some options and run the simuution with the animation in the window. <br>
### Docker
Alternatively use the included Dockerfile to build the docker image. Be aware that running GUI application using docker can be fiddly. To get this running on my mac I followed the instruction [here](https://gist.github.com/paul-krohn/e45f96181b1cf5e536325d1bdee6c949).
