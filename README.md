# Differential-Equations
For usage of given program few things are needed:
- Setup python virtual environment:

   - Enter working directory where this repository was cloned
   - For Linux in command line type:
     ```
     python3 -m venv venv
     ```
      (the last venv can be any name you want your environment to be called)
  -  For Windows
     ```
     python -m venv venv
     ```
- Activate vitrual environment
   - For Linux in comand line type:
     ```
     source venv/bin/activate
     ```
  -  For Windows
     ```
     venv/Scripts/activate
     ```
- Intall package with
     ```
     pip install -e .
     ```
# Considered problem - Heat transport equation
$$-\frac{d}{dx} \left( k(x) \frac{du(x)}{dx} \right) = 100x$$

$$u(2) = 0 ,\hspace{0.5cm}\frac{du(0)}{dx} + u(0) = 20$$

$$
k(x) =
\begin{cases} 
1 & \text{dla } x \in [0, 1] \\
2x & \text{dla } x \in (1, 2]
\end{cases}
$$

Where $$u$$ is unknown function to be found.

$$[0, 2] \ni x \mapsto u(x) \in \mathbb{R}$$

## Solving
Calculation up to the point of determining bilinear form B and linear form L can be found in pdf [here](https://github.com/NCCMNT/Differential-Equations/blob/main/heat%20trasfer%20equation.pdf).
From this point the problem is solved numerically using Finite Element Method in [here](https://github.com/NCCMNT/Differential-Equations/blob/main/main.py).
