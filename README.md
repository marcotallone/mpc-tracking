<a name="readme-top"></a>

<!-- PROJECT SHIELDS -->
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]
[![Gmail][gmail-shield]][gmail-url]

<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href="https://github.com/marcotallone/mpc-tracking">
    <img src="./images/track.png" alt="Logo" width="150" height="150">
  </a>

<h2 align="center">Non-linear Reference Tracking<br>
via Model Predictive Control (MPC)<br>
and Extended Kalman Filter (EKF)</h2>
<h4 align="center">Modelling and Control of Cyber-Physical Systems II Course Exam Project</h4>
<h4 align="center">SDIC Master Degree, University of Trieste (UniTS)</h4>
<h4 align="center">2024-2025</h4>

  <p align="center">
    Implementation of a non-linear Model Predictive Control (MPC) algorithm for reference traking<br>
    with the addition of an Extended Kalman Filter (EKF) for state estimation and robustness to noise.
    <br />
    <br />
    <table>
      <tr>
        <td><a href="//TODO:add link"><strong>Presentation</strong></a></td>
        <td><a href="./report/main.pdf"><strong>Report</strong></a></td>
        <td><a href="./scripts/example.m"><strong>Usage Example</strong></a></td>
      </tr>
    </table>
</div>

<!-- TABLE OF CONTENTS -->
<div align="center">
  <table>
      <tr><td style="text-align: left;">
        <h2>Table of Contents&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</h2>
        <div style="display: inline-block; text-align: left;" align="left">
          <p>
            &nbsp;1. <a href="#author-info">Author Info</a><br>
            &nbsp;2. <a href="#about-the-project">About The Project</a><br>
              &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- <a href="#quick-overview">Quick Overview</a><br>
              &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- <a href="#built-with">Built With</a><br>
              &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- <a href="#project-structure">Project Structure</a><br>
            &nbsp;3. <a href="#getting-started">Getting Started</a><br>
            &nbsp;4. <a href="#usage-examples">Usage Examples</a><br>
            &nbsp;5. <a href="#implementation-details">Implementation Details</a><br>
              &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- <a href="#dynamicalsystem-class">DynamicalSystem Class</a><br>
              &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- <a href="#unicycle-class">Unicycle Class</a><br>
              &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- <a href="#helicopter-class">Helicopter Class</a><br>
              &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- <a href="#mpc-class">MPC Class</a><br>
            &nbsp;6. <a href="#contributing">Contributing</a><br>
            &nbsp;7. <a href="#license">License</a><br>
            &nbsp;8. <a href="#references">References</a><br>
            &nbsp;9. <a href="#acknowledgments">Acknowledgments</a><br>
          </p>
        </div>
      </td></tr>
  </table>
</div>

<!-- AUTHOR INFO-->
## Author Info

| Name | Surname | Student ID | UniTS mail | Google mail | Master |
|:---:|:---:|:---:|:---:|:---:|:---:|
| Marco | Tallone | SM3600002 | <marco.tallone@studenti.units.it> | <marcotallone85@gmail.com> | **SDIC** |

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- ABOUT THE PROJECT -->
## About The Project

>[!WARNING] 
>**Generative Tools Notice**:\
> Generative AI tools have been used as a support for the development of this project. In particular, the [Copilot](https://en.wikipedia.org/wiki/Microsoft_Copilot) generative tool based on [OpenAI GPT 4o](https://en.wikipedia.org/wiki/GPT-4) model has been used as assistance medium in performing the following tasks:
> - writing documentation and comments in the implemented models for `MATLAB`'s `help` function
> - fixing implementation bugs in the `MPC.dense_formulation()` method used to define the QP problem matrices in dense MPC formulation
> - improving variable naming and overall code readability
> - grammar and spelling check both in this README and in the [report](./report/main.pdf)
> - aesthetic improvements in [report](./report/main.pdf) plots
>
> Nevertheless, the auhtor assumes full responsibility for the final content and correctness of the project.

### Quick Overview

Model Predictive Control (MPC) is a control strategy that has been widely
adopted for the control of dynamical systems. The MPC approach is based on the
solution of an optimization problem over a finite horizon, which allows for the
consideration of constraints on the system states and inputs. When used for
reference tracking, the MPC algorithm
computes the optimal control input by minimizing a cost function that penalizes
the deviation of the system states from the desired trajectory.\
However, real-world systems are often affected by noise and disturbances, which
can lead to undesired results with the adoption of an MPC controller. Moreover,
whenever the complete state of a system cannot be fully measured, the use of
state estimators is required to infer the unmeasured components and
apply the MPC algorithm.
To address
this issues, the Extended Kalman Filter (EKF) [<a href="#ref3">3</a>,<a href="#ref4">4</a>]
can be used to estimate the states
of a non-linear system by fusing the available measurements with the system
dynamics and incorporating information about the process and measurement
noise.\
This projects presents the implementation of a Model Predictive Control algorithm
for non-linear dynamical systems based on successive linearization and
discretization of the dynamics around operating points, as well as the
development of an Extended Kalman Filter to estimate the states of the system in
the presence of noise.\
The objectives of this project are both to study the effectiveness of the
non-linear MPC algorithm applied for tracking reference trajectories under
different circumstances and also to reproduce and compare the obtained results
with the ones presented by *Kunz, Huck and Summers* [<a href="#ref5">1</a>]
their work on the tracking of a helicopter model.\
The project is developed in `MATLAB` version `24.2.0 (R2024b)` [<a href="#ref5">5</a>].

<!-- //TODO: add gifs -->
<!-- <div style="text-align: center;">
  <img src="./report/images/visual-words.jpg" alt="Image" width="400"/>
</div> -->

### Project Structure

A complete theoretical description of the implemented models and algorithms can be found formalized in the [report](./report/main.pdf) or pesented in the animated [presentation](//TODO:...) of this repository. Technical implementation details are also provided in section [Implementation Details](#implementation-details) below or can be accessed at any moment using the `help` function in the MATLAB command window.\
The rest of the project is structured as follows:

```bash
.
├── 🖼️ images                       # Project images
├── ⚜️ LICENSE                      # License file
├── 📁 presentation                 # Presentation folder
├── 📜 README.md                    # This file
├── 📁 report                       # Report folder
│   ├── images
│   ├── main.pdf
│   └── ...
├── 📁 results                      # Results folder
│   ├── helicopter_results.csv
│   ├── unicycle_results.csv
│   └── ...
├── 📁 scripts                      # Main MATLAB scripts
│   ├── helicopter_MPC.m
│   ├── helicopter_simulate.m
│   ├── unicycle_MPC.m
│   └── unicycle_simulate.m
└── 📁 src                          # Models implementations
    ├── DynamicalSystem.m
    ├── Helicopter.m
    ├── MPC.m
    └── Unicycle.m
```
  
### Built With

![MATLAB](https://img.shields.io/badge/MATLAB-0076A8?style=for-the-badge&logo=mathworks&logoColor=white)

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- GETTING STARTED -->
## Getting Started

### Requirements

The project is developed in [`MATLAB`](https://www.mathworks.com) version `24.2.0 (R2024b)` [<a href="#ref5">5</a>] and requires the following toolboxes:

- [`Control System Toolbox`](https://www.mathworks.com/products/control.html)
- [`Optimization Toolbox`](https://www.mathworks.com/products/optimization.html)
- [`Symbolic Math Toolbox`](https://www.mathworks.com/products/symbolic.html)
- [`Statistics and Machine Learning Toolbox`](https://www.mathworks.com/products/statistics.html)
- [`Robotics System Toolbox`](https://www.mathworks.com/products/robotics.html)

### Installation

The implemented model classes as well as the MPC class can be found in the [`src/`](./src) folder. In order to run the provided scripts or use the implemented classes in your own scripts, you can simply add the [`src/`](./src) folder to your MATLAB path. This can be done by running the following command in the MATLAB command window:

```matlab
addpath(genpath('path-to-src'))
```

where `'path-to-src'` is the path to the [`src/`](./src) folder relative to your current working directory.\
Alternatively you can also permanently add the [`src/`](./src) folder to your MATLAB path with the following steps:

1. Run 'pathtool' in the MATLAB command window
2. Click on *"Add Folder"* (or *"Add with Subfolders"*) and select the [`src/`](./src) directory.
3. Click *"Save"* to save the path for future MATLAB sessions.

After these steps the [`Unicycle`](./src/Unicycle.m), [`Helicopter`](./src/Helicopter.m) and [`MPC`](./src/MPC.m) classes can be correctly used in any script.

<!-- USAGE EXAMPLES -->
## Usage Examples

The MATLAB scripts in the [`scripts/`](./scripts) folder contain some usage example of the implemented algorithms. In particular the  [`example.m`](./scripts/example.m) script is a step by step tutorial on how to use the implemented classes and methods for the tracking problem. Additionally, the following files can be found with the relative purpose:

- [`unicycle_simulate.m`](./scripts/unicycle_simulate.m): a script to simulate the unicycle model and visualize its time evolution given a determined input sequence
- [`unicycle_MPC.m`](./scripts/unicycle_MPC.m): a script to simulate the unicycle model with the implemented MPC algorithm and the posibility to add process and measurement noise
- [`helicopter_simulate.m`](./scripts/helicopter_simulate.m): a script to simulate the helicopter model and visualize its time evolution given a determined input sequence
- [`helicopter_MPC.m`](./scripts/helicopter_MPC.m): a script to simulate the helicopter model with the implemented MPC algorithm and the posibility to add process and measurement noise

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- IMPLEMENTATION DETAILS -->
## Implementation Details

This section presents some (*not all...*) of the technical details of the implemented classes and methods while theoretical information can be found in the [report](./report/main.pdf).
However, all the implemented classes and methods are also documented respecting `MATLAB`'s documentation conventions. Therefore, further information about any of the implemented models can be found using the `help` function in the MATLAB command window. For example, to get information about the `MPC` class, you can run the following command:

```matlab
help MPC
```

Alternatively, you can also use the `doc` function to open the documentation in the MATLAB Help Browser:

```matlab
doc MPC
```

### `DynamicalSystem` Class

The [`DynamicalSystem`](./src/DynamicalSystem.m) class is an abstract class that defines the interface for all the dynamical systems that can be used with the implemented MPC algorithm. Hence, all the systems that can be implemented have in common the same properties and methods defined in this class. In particular, all systems are defined by the properties:

- `n`: the number of states of the system
- `m`: the number of inputs of the system
- `p`: the number of outputs of the system

And they all share the same implemntation of the `discretize()` method, that uses first order Euler discretization to obtain the discrete-time representation from the linearized continuous-time dynamics of a system around an operating point.\
Then, each class derived from this one must implement the following system-specific methods:

- `dynamics()`: to define the non-linear dynamics of the system
- `simulate()`: to simulate (*possibly using [`ode45`](https://www.mathworks.com/help/matlab/ref/ode45.html)*) the non-linear system dynamics given an input sequence and an initial state
- `output()`: to define the output of the system given the states
- `linearize()`: to linearize the non-linear dynamics around an operating point
- `fix_angles()`: to fix values of reference states angles so that they are defined in the same angular interval (either $[-\pi,+\pi]$ or $[0,2\pi]$) as the measured system states angles (and hence allow the MPC algorithm to compute the correct angular distance)
- `EFK_estimate()`: to estimate the states of the system using the Extended Kalman Filter (EKF)
- `generate_trajectory()`: to generate a reference trajectory for the system to track

### `Unicycle` Class

The [`Unicycle`](./src/Unicycle.m) class is a concrete implementation of the `DynamicalSystem` class that defines the dynamics of a simple unicycle model. The unicycle model is a non-linear system with three states and two inputs, that represents a simple wheeled robot that can move in the plane.\
The model can be initialized by providing the necessary input arguments as in the example below:

```matlab
r  = 0.03;                          % wheel radius
L  = 0.3;                           % distance between wheels
Ts = 0.1;                           % sampling time
x_constraints = [
    -2 2;                           % x position constraints
    -2 2;                           % y position constraints
    -pi 3*pi;                       % heading angle constraints
];
u_constraints = [-50, 50];          % angular velocity constraints
states = 3;                         % number of states
outputs = 2;                        % number of outputs
Q_tilde = 0.75*1e-3*eye(states);    % Process noise covariance
R_tilde = 1e-2*eye(outputs);        % Measurement noise covariance
P0 = eye(states);                   % Initial state covariance

model = Unicycle(r, L, Ts, x_constraints, u_constraints, P0, Q_tilde, R_tilde);
```

### `Helicopter` Class

The [`Helicopter`](./src/Helicopter.m) class is a concrete implementation of the `DynamicalSystem` class that defines the dynamics of a simple helicopter model. The helicopter model is a non-linear system with $8$ states and $4$ inputs, that represents a simple helicopter that can move in the three-dimensional space.\
The model can be initialized by providing the necessary input arguments as in the example below:

```matlab
bx = 2; by = 2; bz = 18; bpsi = 111;
kx = -0.5; ky = -0.5; kpsi = -5; ki = 2;
parameters = [bx; by; bz; bpsi; kx; ky; kpsi; ki];
g = 9.81;
Ts = 0.1;                           % sampling time
x_constraints = [                   % state constraints
    -10, 10;                        % xI constraints
    -10, 10;                        % yI constraints
    -10, 10;                        % zI constraints
    -10, 10;                        % vxB constraints
    -10, 10;                        % vyB constraints
    -10, 10;                        % vzB constraints
    -pi, 3*pi;                      % psi constraints
    -25, 25;                        % omega constraints
];
u_constraints = [                   % input constraints
    -1, 1;                          % ux constraints
    -1, 1;                          % uy constraints
    -1, 1;                          % uz constraints
    -1, 1;                          % upsi constraints            
];
states = 8;                         % number of states
outputs = 4;                        % number of outputs
Q_tilde = 0.75*1e-3*eye(states);    % Process noise covariance
R_tilde = 1e-2*eye(outputs);        % Measurement noise covariance
P0 = eye(states);                   % Initial state covariance

model = Helicopter(parameters, Ts, x_constraints, u_constraints, P0, Q_tilde, R_tilde);
```

### `MPC` Class

The [`MPC`](./src/MPC.m) class defines the Model Predictive Control algorithm for non-linear dynamical systems. The MPC algorithm is based on the successive linearization and discretization of the dynamics around an operating point, and the solution of a Quadratic Programming (QP) problem to compute the optimal control input at each discretized time step.\
The `MPC` class object can be initialized providing the fllowing input arguments depending on the conditions on which the MPC algorithm is applied:

```matlab
mpc = MPC(model, x0, Tend, N, Q, R, x_ref, u_ref, preview, formulation, noise, debug)
```

where:

- `model`: the dynamical system model to control that inherits from the `DynamicalSystem` class
- `x0`: the initial state of the system
- `Tend`: the final time of the simulation
- `N`: the prediction horizon of the MPC algorithm
- `Q`: the state cost matrix of the MPC algorithm
- `R`: the input cost matrix of the MPC algorithm
- `x_ref`: the reference trajectory to track
- `u_ref`: the reference input sequence to track
- `preview`: boolean flag: if enabled (i.e. `1`, default choice) the MPC algorithm sees N steps ahead in the reference trajectory, otherwise (`0`) it sees only the current step reference N times
- `formulation`: the formulation of the QP problem to solve: either `'0: dense'` or `'1: sparse'`
- `noise`: boolean flag: if enabled (i.e. `1`) the MPC algorithm considers process and measurement noise in the system and the `EFK_estimate()` method is called at each time step, otherwise (`0`, default choice) the MPC algorithm considers a noise-free system
- `debug`: boolean flag: if enabled (i.e. `1`) the MPC algorithm prints debug information at each time step, otherwise (`0`, default choice) it does not print any information

>[!NOTE]
> The `MPC` will internally look at the discretization time step `Ts` of the dynamical system model to compute the total number of iterations of the MPC algorithm. Precisely, the total number `Nsteps` of iterations is computed as
> ```matlab
>t = 0:model.Ts:Tend;
>Nsteps = length(t)-(N+1);
>```

After initialization with the desired setup, one can then collect the states of the dynamical system during simulation as well as the optimal input sequence at each time step by running the `optimize()` method of the `MPC` object.

```matlab
[x, u] = mpc.optimize();
```

At each time step, the `optimize()` method internally uses the [`quadprog`](https://www.mathworks.com/help/optim/ug/quadprog.html) function to solve the QP problem resultig from either the *dense* or *sparse* formulation chosen during initialization.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- CONTRIBUTING -->
## Contributing

The goal of this repository was to implement a non-linear MPC algorithm with the addition of the Extended Kalman Filter state estimator and to reproduce the results presented in the referenced papers, in the context of a university exam project. However, if you have a suggestion that would make this better or extend its functionalities and want to share it with me, please fork the repo and create a pull request. You can also simply open an issue with the tag "enhancement" or "extension".\
Suggested contribution procedure:

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- LICENSE -->
## License

Distributed under the MIT License. See [`LICENSE`](./LICENSE) for more information.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- REFERENCES -->
## References

<a id="ref1"></a>
[1] K. Kunz, S. M. Huck, T. H. Summers, *"Fast Model Predictive Control of miniature helicopters"*, 2013 European Control Conference (ECC), 2013, Pages 1377-1382, [https://doi.org/10.23919/ECC.2013.6669699](https://ieeexplore.ieee.org/document/6669699)

<a id="ref2"></a>
[2] R.M. Murray, *"Optimization Based Control"*, California Institute of Technology, 2023, Chapter Trajectory Generation and Differential Flatness, [http://www.cds.caltech.edu/~murray/books/AM08/pdf/obc-complete_12Mar2023.pdf](http://www.cds.caltech.edu/~murray/books/AM08/pdf/obc-complete_12Mar2023.pdf)

<a id="ref3"></a>
[3] R. E. Kalman, *"A New Approach to Linear Filtering and Prediction Problems"*, Journal of Basic Engineering, Volume 82, 1960, Pages 35-45, [https://doi.org/10.1115/1.3662552](https://asmedigitalcollection.asme.org/fluidsengineering/article-pdf/82/1/35/5518977/35_1.pdf)

<a id="ref4"></a>
[4] R. E. Kalman, R. S. Bucy, *"New Results in Linear Filtering and Prediction Theory"*, Journal of Basic Engineering, Volume 83, 1961, Pages 95-108, [https://doi.org/10.1115/1.3658902](https://asmedigitalcollection.asme.org/fluidsengineering/article-pdf/83/1/95/5503549/95_1.pdf)

<a id="ref5"></a>
[5] The MathWorks Inc., *"MATLAB version: 24.2.0 (R2024b)"*, The MathWorks Inc., 2024, [https://www.mathworks.com](https://www.mathworks.com)

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- ACKNOWLEDGMENTS -->
## Acknowledgments

- [Modelling and Control of Cyberphysical Systems II Course (UniTS, Fall 2024)](https://dsai.units.it/index.php/courses/modelling-and-control-of-cyber-physical-systems-ii/) (*access restricted to UniTS students and staff*)
- [Best-README-Template](https://github.com/othneildrew/Best-README-Template?tab=readme-ov-file): for the README template
- [Flaticon](https://www.flaticon.com/free-icon/track_7171906): for the icons used in the README

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[forks-shield]: https://img.shields.io/github/forks/marcotallone/mpc-tracking.svg?style=for-the-badge
[forks-url]: https://github.com/marcotallone/mpc-tracking/network/members
[stars-shield]: https://img.shields.io/github/stars/marcotallone/mpc-tracking.svg?style=for-the-badge
[stars-url]: https://github.com/marcotallone/mpc-tracking/stargazers
[issues-shield]: https://img.shields.io/github/issues/marcotallone/mpc-tracking.svg?style=for-the-badge
[issues-url]: https://github.com/marcotallone/mpc-tracking/issues
[license-shield]: https://img.shields.io/github/license/marcotallone/mpc-tracking.svg?style=for-the-badge
[license-url]: https://github.com/marcotallone/mpc-tracking/blob/master/LICENSE
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-blue?style=for-the-badge&logo=linkedin&logoColor=white&colorB=0077B5
[linkedin-url]: https://linkedin.com/in/marco-tallone-40312425b
[gmail-shield]: https://img.shields.io/badge/-Gmail-red?style=for-the-badge&logo=gmail&logoColor=white&colorB=red
[gmail-url]: mailto:marcotallone85@gmail.com
