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

Incident cosmic ray muons were simulated in a monte carlo fashion, predominantly making use of the random function from the Numpy random library. Here we are interested in mapping the uniform probability distribution on the interval \[0, 1) provided by the random function to a general pdf P(y) on some specified interval between y<sub>0</sub> and y<sub>1</sub>. In order for this to work we require that the CDFs (cumulative density functions) are equal. This means that for some randomly generated xin from the random function we have the following equality <br>
![equation](https://latex.codecogs.com/gif.latex?%5Cint_%7Bx_0%7D%5E%7Bx_%7Bin%7D%7D%20u%28x%29%20dx%20%3D%20%5Cint_%7By_0%7D%5E%7By_%7Bout%7D%7D%20P%28y%29%20dy)
 <br> where u(x) is the uniform distribution provided by the random function and yout is the desired output. The LHS of this comes out as just xin which leads to the following equation for yout <br>
![equation](https://latex.codecogs.com/gif.latex?y_%7Bout%7D%20%3D%20F%5E%7B-1%7D%28x_%7Bin%7D%29)
<br>
## Numerical simulation of ion drift diffusion 
