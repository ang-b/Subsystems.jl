# Subsystems

A low-expectation package to work with N interconnected systems of the form 

<img src="https://render.githubusercontent.com/render/math?math=%5Cdot%7Bx%7D_i%20%3D%20f(x_i%2C%20x_j%2C%20u_i)">

for i,j = 1, ..., N.
The code only implements discrete-time equations and it is very experimental and unstable as I use it for my own research and applications. 
The package also includes observer implementations for classic Unknown Input Observers [1], Unbiased Minimum Variance Estimators [2,3], and Distributed Luenberger Observers [4].

### Disclaimer
I am a Julia newcomer and the code herein is the result of rushed fiddling to get simulations done as quickly as possible. 
I do intend to improve and maintain this because it is lots of fun, so any constructive criticism is welcome.

### References 

[1] Chen, J., Patton, R. J., & Zhang, H. Y. (1996). Design of unknown input observers and robust fault detection filters. International Journal of control, 63(1), 85-10

[2] Gillijns, S., & De Moor, B. (2007). Unbiased minimum-variance input and state estimation for linear discrete-time systems. Automatica, 43(1), 111-116

[3] Barboni, A., Gallo, A. J. et al. (2019). A Distributed Approach for the Detection of Covert Attacks in Interconnected Systems with Stochastic Uncertainties. In 2019 IEEE 58th Conference on Decision and Control (CDC). IEEE.

[4] Barboni, A., Rezaee, H. et al. (2020). Detection of covert cyber-attacks in interconnected systems: a distributed model-based approach. IEEE Transactions on Automatic Control, 65(9), 3728-3741.
