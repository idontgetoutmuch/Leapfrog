% Planetary Simulation with Excursions in Symplectic Manifolds
% Dominic Steinitz
% 7th June 2013

This article attempts to show that Haskell [@Marlow_haskell2010]
performs reasonably well on numerical problems.

When I started to do this, it seemed straightforward enough: pick a
problem which admitted a numerical solution, find an algorithm and
code it up. I chose the problem of orbital dynamics as I had always
been fascinated by the precession of the perihelion of Mercury (which
is mainly caused by the pull of the other planets in the Solar System)
and because this admits of at least two different methods of numerical
solution both of which I hope will show the power of Haskell in this
area. This led to the selection of algorithm and I read that one
should prefer a symplectic method such as the Leapfrog which conserves
the energy of a system (a highly desirable requirement when modelling
orbital dynamics). My conscience would not let me write about such a
method without being able to explain it. This led into the Hamiltonian
formulation of classical mechanics, symplectic manifolds and
symplectic (numerical) methods.

The reader interested in the Haskell implementations and performance
comparisons with other programming languages can read the introduction
and skip to ????. I apologise in advance to experts in classical
mechanics, symplectic geometery and numerical analysis and can only
hope I have not traduced their subjects too much.

Introduction
------------

Forget about Newton and suppose you are told that the way to do
mechanics is to write down the total energy of the system in which you
are interested and then apply Hamilton's equations.

Consider a mass of $m$ attached to a light rod of length $l$ which is
attached to a point from which it can swing freely in a plane. Then
the kinetic energy is:

$$
\frac{1}{2}mv^2 = \frac{1}{2}ml^2\dot{\theta}^2
$$

and the potential energy (taking this to be 0 at $\theta = 0$) is:

$$
mgl(1 - \cos\theta)
$$

Thus the Hamiltonian is:

$$
\cal{H} = \frac{1}{2}ml^2\dot{\theta}^2 + mgl(1 - \cos\theta)
$$

Let us set the generalized momentum

$$
p = \frac{\partial\cal{L}}{\partial\dot{\theta}} = ml^2\dot{\theta}
$$

Then we can re-write the Hamiltonian as:

$$
\cal{H} = \frac{p^2}{2ml^2} + mgl(1 - \cos\theta)
$$

Applying Hamilton's equations we obtain

$$
\begin{aligned}
\dot{\theta} &=  \frac{\partial\cal{H}}{\partial p}      = \frac{p}{ml^2} \\
\dot{p}      &= -\frac{\partial\cal{H}}{\partial \theta} = -mgl\sin\theta
\end{aligned}
$$

Differentiating the first equation with respect to time we then obtain
the familiar equation describing the motion of a simple pendulum.

$$
\ddot{\theta} = \frac{\dot{p}}{ml^2} = \frac{-mgl\sin\theta}{ml^2} = -\frac{g}{l}\sin\theta
$$

Now we would like to calculate the pendulum's position and velocity at
a given time. The obvious starting place is to use the explicit Euler
method.

$$
\begin{aligned}
\theta_{n+1} &=  h\frac{p_n}{ml^2} \\
p_{n+1}      &= -hmgl\sin\theta_n
\end{aligned}
$$

{-# OPTIONS_GHC -Wall                     #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing  #-}
{-# OPTIONS_GHC -fno-warn-type-defaults   #-}

> {-# LANGUAGE TupleSections                #-}
> {-# LANGUAGE NoMonomorphismRestriction    #-}

> module Symplectic (
>     blsEE
>   , tls
>   , bls
>   , trs
>   , brs
>   ) where

> stepMomentumEE :: Double -> Double -> Double -> Double -> Double
> stepMomentumEE m l p q = p -  h * m * g * l * sin q

> stepPositionEE :: Double -> Double -> Double -> Double -> Double
> stepPositionEE m l p q = q + h * p / (m * l^2)

> stepOnceEE :: Double -> Double -> Double -> Double -> (Double, Double)
> stepOnceEE m l p q = (newP, newQ)
>   where
>     newP = stepMomentumEE m l p q
>     newQ = stepPositionEE m l p q

> h, m, l, g :: Double
> h = 0.01  -- Seconds
> l = 1.0  -- Metres
> m = 1.0  -- Kilograms
> g = 9.81 -- Metres * Seconds^-2

> initTheta, initThetaDot, initP :: Double
> initTheta    = 0.0
> initThetaDot = 1.7
> initP        = m * l^2 * initThetaDot

> runEE :: Double -> Double -> [(Double, Double)]
> runEE initP initTheta = iterate (uncurry (stepOnceEE m l)) (initP, initTheta)

```{.dia width='800'}
import Symplectic
import SymplecticDia

diaEE :: DiagramC
diaEE = test tickSize [ (cellColourEE0, take nPlotPoints $ blsEE)
                      ]
dia = diaEE
```

As we can see from the diagram above, energy is not conserved but
increases steadily over time, an undesirable state of affairs.

Instead let us apply the the symplectic Euler method:

$$
\begin{aligned}
p_{n+1} = p_n - hmgl\sin\theta_n \\
\theta_{n+1} = \theta_n + \frac{hp_{n+1}}{2ml^2}
\end{aligned}
$$

> stepMomentum :: Double -> Double -> Double -> Double -> Double
> stepMomentum m l p q = p -  h * m * g * l * sin q

> stepPosition :: Double -> Double -> Double -> Double -> Double
> stepPosition m l p q = q + h * p / (m * l^2)

> stepOnce :: Double -> Double -> Double -> Double -> (Double, Double)
> stepOnce m l p q = (newP, newQ)
>   where
>     newP = stepMomentum m l p q
>     newQ = stepPosition m l newP q


```{.dia width='800'}
import Symplectic
import SymplecticDia

dia' :: DiagramC
dia' = test tickSize [ (cellColour0, take nPlotPoints $ bls)
                     ]

dia = dia'
```

In this case the energy is conserved so this looks like a good
candidate for simulating orbital dynamics. But why does this work? It really looks very similar to the explicit Euler method.

Theory
------

```{.dia width='800'}
illustrateBezier c0 c1 c2 c3 x2 x3
    = endpt  # translate x3
    <> l3a
    <> fromSegments [bezier3 c3 c0 x3, bezier3 c1 c2 x2]
  where
    dashed  = dashing [0.1,0.1] 0
    endpt   = circle 0.05 # fc red  # lw 0
    l3a     = fromOffsets [r2 (1, 2)] # translate (r2 (1, 1)) # dashed

x2      = r2 (3,-1) :: R2         -- endpoint
x3      = r2 (1, 1) :: R2         -- endpoint
[c0,c1,c2,c3,c4] = map r2 [(-1, -3), (1,2), (2,0), (-3,0), (-1, -2)]   -- control points

example = illustrateBezier c0 c1 c2 (-c2) x2 x3
dia = example
```

$$
\begin{aligned}
p_{n+1} &= p_n - h\nabla_q H(p_{n+1}, q_n) \\
q_{n+1} &= q_n + h\nabla_p H(p_{n+1}, q_n)
\end{aligned}
$$

We check that this really is symplectic. First suppose we have two functions:

$$
\begin{aligned}
x &= u - f(x,v) \\
y &= v + g(x,v) \\
\end{aligned}
$$

Then we can find partial derivatives:

$$
\begin{aligned}
dx &= du - \frac{\partial f}{\partial x}dx - \frac{\partial f}{\partial v}dv \\
dy &= dv + \frac{\partial g}{\partial x}dx + \frac{\partial g}{\partial v} dv \\
\end{aligned}
$$

$$
\begin{aligned}
\frac{\partial x}{\partial u} &= 1 - \frac{\partial f}{\partial x}\frac{\partial x}{\partial u} \\
\frac{\partial x}{\partial v} &= -\frac{\partial f}{\partial x}\frac{\partial x}{\partial v} -\frac{\partial f}{\partial v} \\
\frac{\partial y}{\partial u} &= \frac{\partial g}{\partial x}\frac{\partial x}{\partial u} \\
\frac{\partial y}{\partial v} &= 1 + \frac{\partial g}{\partial x}\frac{\partial x}{\partial v} + \frac{\partial g}{\partial v}
\end{aligned}
$$

Re-arranging:

$$
\begin{aligned}
\frac{\partial x}{\partial u}(1 + \frac{\partial f}{\partial x}) &= 1 \\
\frac{\partial x}{\partial v}(1 + \frac{\partial f}{\partial x}) &= -\frac{\partial f}{\partial v} \\
\frac{\partial y}{\partial u} -\frac{\partial g}{\partial x}\frac{\partial x}{\partial u} &= 0 \\
\frac{\partial y}{\partial v} - \frac{\partial g}{\partial x}\frac{\partial x}{\partial v} &= 1 + \frac{\partial g}{\partial v}
\end{aligned}
$$

Pulling everything together in matrix form:

$$
\begin{bmatrix}
1 + \frac{\partial f}{\partial x} & 0 \\
-\frac{\partial g}{\partial x} & 1
\end{bmatrix}
\,
\begin{bmatrix}
\frac{\partial x}{\partial u} & \frac{\partial x}{\partial v} \\
\frac{\partial y}{\partial u} & \frac{\partial y}{\partial v}
\end{bmatrix}
=
\begin{bmatrix}
1 & -\frac{\partial f}{\partial v} \\
0 & 1 + \frac{\partial g}{\partial v}
\end{bmatrix}
$$



> runSE :: Double -> Double -> [(Double, Double)]
> runSE initP initTheta = iterate (uncurry (stepOnce m l)) (initP, initTheta)

> bls   = runSE initP         initTheta
> blsEE = runEE initP         initTheta
> brs   = runSE (initP + 1.0) initTheta
> trs   = runSE (initP + 1.0) (initTheta + 1.0)
> tls   = runSE initP         (initTheta + 1.0)
>
> areaParGram (x1, y1) (x2, y2) = (x2 - x1) * (y2 - y1)
> areas = zipWith areaParGram bls trs
>
> hamiltonian :: Double -> Double -> Double -> Double -> Double
> hamiltonian m l p q = (p^2 / (2 * m * l^2)) + (m * g * l * (1 - cos q))


```{.dia width='800'}
import Symplectic
import SymplecticDia

dia' :: DiagramC
dia' = test tickSize [ (cellColour0, take nPlotPoints $ bls)
                     , (cellColourEE0, take nPlotPoints $ blsEE)
                     , (cellColour1, take nPlotPoints $ brs)
                     , (cellColour2, take nPlotPoints $ trs)
                     , (cellColour3, take nPlotPoints $ tls)
                     ]

dia = dia'
```

Performance
-----------

Bibliography
------------

