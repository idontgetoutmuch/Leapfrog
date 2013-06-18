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

> {-# LANGUAGE NoMonomorphismRestriction    #-}
> {-# LANGUAGE FlexibleContexts             #-}
> {-# LANGUAGE ScopedTypeVariables          #-}

> module Symplectic (
>     blsEE
>   , tls
>   , bls
>   , trs
>   , brs
>   , simPlanets
>   ) where

> import Data.Array.Repa hiding ((++), zipWith)
> import qualified Data.Array.Repa as Repa
> import Data.Array.Repa.Algorithms.Matrix
> import Control.Monad
> import Control.Monad.Identity
> import Text.Printf

> import Debug.Trace

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

The Canonical Symplectic Form for the Cotangent Bundle
======================================================

The cotangent bundle has a canonical symplectic 2-form and hence is a symplectic manifold.

Let $\pi : T^* M \longrightarrow M$ be the projection function from
the cotangent bundle to the base manifold, that is, $\pi(x,\xi) =
x$. Then $\pi_* : T(T^*M) \longrightarrow TM$ and we can define a 1-form (the
canonical or tautological 1-form) on $v \in T_{(x,\xi)}(T^* M)$ as

$$
\theta_{(x,\xi)} (v) = \xi(\pi_* v)
$$

By definition:

$$
\pi_* \bigg(\frac{\partial}{\partial x_i}\bigg)(f) = \frac{\partial}{\partial x_i}(f \circ \pi) = \frac{\partial f}{\partial x_i}
$$

and

$$
\pi_* \bigg(\frac{\partial}{\partial \xi_i}\bigg)(f) = \frac{\partial}{\partial \xi_i}(f \circ \pi) = 0
$$

If we then write $v \in T_{(x,\xi)}(T^*M)$ in co-ordinate form:

$$
v = a^i\frac{\partial}{\partial x_i} + \alpha^i\frac{\partial}{\partial \xi_i}
$$

we have

$$
\pi_*v = a^i\pi_*\frac{\partial}{\partial x_i} + \alpha^i\pi_*\frac{\partial}{\partial \xi_i} = a^i\frac{\partial}{\partial x_i}
$$

Taking $\xi = \xi_i dx^i$ we have that

$$
\xi(\pi_*v) = \xi_i a^i
$$

Thus we have:

$$
\theta_{(x,\xi)} = \xi_i dx^i
$$

We then have a closed 2-form

$$
\omega = d\theta
$$

In co-ordinate terms:

$$
\omega = d\theta = d (\xi_i dx^i) = d\xi^i \wedge dx^i
$$

Hamiltonian Vector Fields
=========================

Without proof let us record the following fact. Let $(M, \omega)$ be a
symplectic manifold. Then there exists a bundle isomorphism
$\tilde{\omega} : TM \longrightarrow T^*M$ defined by
$\tilde{\omega}(X_p)(Y_p) = \omega_p(X_p, Y_p)$.

This is analagous to the isomorphism one can derive in a (semi)
Riemannian manifold with the metric in some sense playing the role of
the 2-form (see [@o1983semi] for example).

Now suppose we have a smooth function $H : M \longrightarrow
\mathbb{R}$ then we can form the 1-form $dH$ and we can define the
Hamiltonian vector field $X_H = \tilde{\omega}^{-1}(dH)$.

In co-ordinates we have:

$$
\omega\bigg(\frac{\partial}{\partial x_i}, \frac{\partial}{\partial x_j}\bigg) = d\xi^i \wedge dx^i\bigg(\frac{\partial}{\partial x_i}, \frac{\partial}{\partial x_j}\bigg) = d\xi^i\bigg(\frac{\partial}{\partial x_i}\bigg)dx^i\bigg(\frac{\partial}{\partial x_j}\bigg) - d\xi^i\bigg(\frac{\partial}{\partial x_j}\bigg)dx^i\bigg(\frac{\partial}{\partial x_i}\bigg) = 0
$$

Euler Symplectic
================

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

Planetary Motion
----------------

Normally we would express the gravitational constant in SI units but
to be consistent with [@hairer2010geometric] we use units in which
distances are expressed in astronomical units, masses are measured
relative to the sun and time is measured in earth days.

$$
{\cal H} = \frac{1}{2}\sum_{i=0}^n \frac{p_i^\top p_i}{m_i} - \frac{G}{2}\sum_{i=0}^n\sum_{j \neq i} \frac{m_i m_j}{\|q_i - q_j\|}
$$

Applying Hamilton's equations we obtain

$$
\begin{aligned}
\dot{q_k^a} &=  \frac{\partial\cal{H}}{\partial p_k^a} = \frac{p_k^a}{m_k} \\
\dot{p_k^a} &= -\frac{\partial\cal{H}}{\partial q_k^a} = G\sum_{j \neq k}m_k m_i \frac{q_k^a - q_j^a}{\|q_k - q_j\|^3}
\end{aligned}
$$

In this case it easy to see that these are the same as Newton's laws of motion.

Applying the Euler symplectic method we obtain:

$$
\begin{aligned}
q_k^{n+1} &= q_k^n + h \frac{p_k^n}{m_k} \\
p_k^{n+1} &= p_k^n + h G\sum_{j \neq k}m_k m_i \frac{q_k^{n+1} - q_j^{n+1}}{\|q_k^{n+1} - q_j^{n+1}\|^3}
\end{aligned}
$$

> stepPositionP :: ( Monad m
>                  , Source a Double
>                  , Source b Double
>                  , Source c Double
>                  ) =>
>                  Double ->
>                  Array a DIM2 Double ->
>                  Array b DIM1 Double ->
>                  Array c DIM2 Double ->
>                  m (Array U DIM2 Double)
> stepPositionP h qs ms ps = do
>   computeP $ qs +^ (ps *^ h2 /^ ms2)
>     where
>       (Z :. i :. j) = extent ps
>
>       h2  = extend (Any :. i :. j) $ fromListUnboxed Z [h]
>       ms2 = extend (Any :. j) ms

> stepMomentumP :: forall a b c m . ( Monad m
>                  , Source a Double
>                  , Source b Double
>                  , Source c Double
>                  ) =>
>                  Double ->
>                  Array a DIM2 Double ->
>                  Array b DIM1 Double ->
>                  Array c DIM2 Double ->
>                  m (Array U DIM2 Double)
> stepMomentumP h qs ms ps =
>   do fs <- sumP $ transpose $ zeroDiags' fss
>      let foo :: Array D DIM2 Double
>          foo = qs +^ qs -^ qs
>      bar <- computeP foo :: m (Array U DIM2 Double)
>      -- trace ("qs:\n" ++ show bar) $ return ()
>      let foo :: Array D DIM2 Double
>          foo = ps +^ ps -^ ps
>      bar <- computeP foo :: m (Array U DIM2 Double)
>      -- trace ("ps:\n" ++ show bar) $ return ()
>      -- trace ("fs:\n" ++ show fs) $ return ()
>      urk <- computeP (fs *^ dt2 /^ ms2) :: m (Array U DIM2 Double)
>      -- trace ("Updates:\n" ++ show urk) $ return ()
>      baz <- computeP (ps +^ (fs *^ dt2 {- /^ ms2 -}))
>      -- trace ("newPs:\n" ++ show baz) $ return baz
>      return baz
>   where
>     is = repDim2to3Outer $ prodPairsMasses ms
>     qDiffs = pointDiffs qs
>     preDs = Repa.map (^3) $
>             Repa.map sqrt $
>             sumS $
>             Repa.map (^2) $
>             qDiffs
>     ds    = repDim2to3Outer preDs
>     preFs = Repa.map (* (negate gConst)) $
>             qDiffs /^ ds
>     fss = is *^ preFs
>     Z :.i :. j :. k = extent fss
>     dt2             = extend (Any :. i :. k) $ fromListUnboxed Z [h]
>     ms2             = extend (Any :. k) ms

> stepOnceP :: ( Monad m
>              , Source a Double
>              , Source b Double
>              , Source c Double
>              ) =>
>              Double ->
>              Array b DIM1 Double ->
>              Array a DIM2 Double ->
>              Array c DIM2 Double ->
>              m (Array U DIM2 Double, Array U DIM2 Double)
> stepOnceP h ms qs ps = do
>   newPs <- stepMomentumP h qs ms ps
>   -- trace ("newPs:\n" ++ show newPs) $ return ()
>   newQs <- stepPositionP h qs ms newPs
>   -- trace ("newQs:\n" ++ show newQs) $ return ()
>   return (newQs, newPs)

> -- FIXME
> repDim2to3Outer a = extend (Any :. spaceDim) a



> gConst :: Double
> -- gConst = 1.0 -- 2.95912208286e-4
> gConst = 6.67384e-11        -- gravitational constant



> spaceDim :: Int
> spaceDim = 3

> type Momenta   = Array U DIM2 Double
> type Positions = Array U DIM2 Double
> type Masses    = Array U DIM1 Double
>
> zeroDiags x = traverse x id f
>   where
>     f _ (Z :. i :. j) | i == j    = 0.0
>                       | otherwise = x!(Z :. i :. j)
>                                     
> zeroDiags' x = traverse x id f
>   where
>     f _ (Z :. i :. j :. k) | i == j    = 0.0
>                            | otherwise = x!(Z :. i :. j :. k)

> hamiltonianP :: Masses -> Momenta -> Positions -> IO Double
> hamiltonianP ms ps qs = do preKes <- sumP $ ps *^ ps
>                            ke     <- sumP $ preKes /^ ms
>
>                            ds2 <- sumP $ Repa.map (^2) $ pointDiffs qs
>                            let ds   = Repa.map sqrt ds2
>                                is   = prodPairsMasses ms
>                                pess = zeroDiags $ Repa.map (* (negate gConst)) $ is /^ ds
>                            pes <- sumP pess
>                            pe  <- sumP pes
>                            te :: Array U DIM0 Double <- computeP $ ke +^ pe
>                            return $ head $ toList $ Repa.map (* 0.5) te

> repDim1To2Outer :: Source a Double =>
>                    Array a DIM1 Double ->
>                    Array D DIM2 Double
> repDim1To2Outer a = extend (Any :. i :. All) a
>   where (Z :. i) = extent a

> prodPairsMasses :: Source a Double =>
>                    Array a DIM1 Double ->
>                    Array D DIM2 Double
> prodPairsMasses ms = ns *^ (transpose ns)
>   where
>     ns = repDim1To2Outer ms

> transposeOuter :: Source a Double =>
>               Array a DIM3 Double ->
>               Array D DIM3 Double
> transposeOuter qs = backpermute (f e) f qs
>   where
>     e = extent qs
>     f (Z :. i :. i' :. j) = Z :. i' :. i :. j
>
> pointDiffs :: Source a Double =>
>               Array a DIM2 Double ->
>               Array D DIM3 Double
> pointDiffs qs = qss -^ (transposeOuter qss)
>   where qss = replicateRows qs
>
> replicateRows :: Source a Double =>
>                  Array a DIM2 Double ->
>                  Array D DIM3 Double
> replicateRows a = extend (Any :. i :. All) a
>   where (Z :. i :. _j) = extent a
>

> k :: Double
> k = 24*60*60                -- seconds in a day

> days, t, timestepDays, dt :: Double
> days = 36500*100            -- total time in days
> t = days*k                  -- total time
> timestepDays = 10           -- timestep in days
> dt = timestepDays*k         -- timestep

> type Distance = Double
> type Mass     = Double
> type Force    = Double
> type Speed    = Double
> type Energy   = Double

Now we need some initial conditions to start our simulation.

  [jupiter]: http://en.wikipedia.org/wiki/Jupiter

> sunMass, jupiterMass, earthMass :: Mass
> sunMass     = 1.9889e30
> jupiterMass = 1.8986e27
> earthMass   = 5.9722e24

> jupiterAphelion   :: Distance
> jupiterAphelion   = 8.165208e11
> jupiterPerihelion :: Distance
> jupiterPerihelion = 7.405736e11
> jupiterEccentrity :: Double     -- Eccentricity is dimensionless
> jupiterEccentrity = 4.8775e-2
>
> jupiterMajRad :: Distance
> jupiterMajRad = (jupiterPerihelion + jupiterAphelion) / 2

> earthAphelion   :: Distance
> earthAphelion   = 1.520982e11
> earthPerihelion :: Distance
> earthPerihelion = 1.470983e11
> earthEccentrity :: Double     -- Eccentricity is dimensionless
> earthEccentrity = 1.6711e-2
>
> earthMajRad :: Distance
> earthMajRad = (earthPerihelion + earthAphelion) / 2

Kepler's third law states, "The square of the orbital period of a
planet is directly proportional to the cube of the semi-major axis of
its orbit". Here it is in mathematical form:

$$
T^2 = \frac{4 \pi^2 a^3}{GM}
$$

where $T$ is the period of the orbit, $a$ is the major radius of the
elliptical orbit (Kepler's first law: "The orbit of every planet is an
ellipse with the Sun at one of the two foci"), $G$ is the
gravitational constant and $M$ is the mass of the sun.

From this we can calculate the mean angular velocity: $n = 2\pi / T$.

> nJupiter :: Double
> nJupiter = sqrt $ gConst * sunMass / jupiterMajRad^3

$$
\begin{align*}
r &= \frac{a(1 - e^2)}{1 - e\cos\theta} \\
r^2\dot{\theta} &= \sqrt{(1 - e^2)}na^2 \\
GM_{\rm Sun} &= n^2a^3
\end{align*}
$$

where $G$ is the gravitational constant, $n = \frac{2\pi}{T}$ is the
mean angular orbital velocity, $a$ is the major access of the planet's
ellipse and $e$ is the eccentricity.

Finally we can calculate Jupiter's velocity by assuming that its
perihelion is on the $x$-axis and that its velocity in the $x$
direction must be $0$.

Let us calculate the initial conditions assuming that Jupiter starts
at its perihelion. The angular velocity at that point is entirely in
the negative $y$ direction.

With the Leapfrog Method we need the velocity to be half a time step
before the perihelion.

If we take $(0.0, -v_p, 0.0)$ to be the velocity of Jupiter at the
perihelion then if $\delta\theta$ is the angle with respect to the
negative y-axis at half a time step before Jupiter reaches the
perihelion then the velocity of Jupiter at this point is given by
simple trigonometry:
$$
(-v_p \sin(\delta\theta), -v_p \cos(\delta\theta), 0.0) \approx (-v_p\delta\theta, -v_p(1-\delta\theta^2 / 2), 0.0)
$$

> jupiterThetaDotP :: Double -- radians per second
> jupiterThetaDotP = nJupiter *
>                    jupiterMajRad^2 *
>                    sqrt (1 - jupiterEccentrity^2) / jupiterPerihelion^2
> jupiterDeltaThetaP :: Double -- radians
> jupiterDeltaThetaP = jupiterThetaDotP * dt / 2
>
> jupiterVPeri :: Speed
> jupiterVPeri = jupiterThetaDotP * jupiterPerihelion
>
> jupiterInitX :: Speed
> jupiterInitX = negate $ jupiterVPeri * jupiterDeltaThetaP
>
> jupiterInitY :: Speed
> jupiterInitY = negate $ jupiterVPeri * (1 - jupiterDeltaThetaP^2 / 2)
>
> jupiterV :: (Speed, Speed, Speed)
> jupiterV = (jupiterInitX, jupiterInitY, 0.0)
>
> jupiterR :: (Distance, Distance, Distance)
> jupiterR = (negate jupiterPerihelion, 0.0, 0.0)

We can do the same for Earth but we assume the earth is at its
perihelion on the opposite side of the Sun to Jupiter.

> nEarth :: Double
> nEarth = sqrt $ gConst * sunMass / earthMajRad^3
>
> earthThetaDotP :: Double -- radians per second
> earthThetaDotP = nEarth *
>                  earthMajRad^2 *
>                  sqrt (1 - earthEccentrity^2) / earthPerihelion^2
> earthDeltaThetaP :: Double -- radians
> earthDeltaThetaP = earthThetaDotP * dt / 2
>
> earthVPeri :: Speed
> earthVPeri = earthThetaDotP * earthPerihelion
>
> earthInitX :: Speed
> earthInitX = earthVPeri * earthDeltaThetaP
>
> earthInitY :: Speed
> earthInitY = earthVPeri * (1 - earthDeltaThetaP^2 / 2)
>
> earthV :: (Speed, Speed, Speed)
> earthV = (earthInitX, earthInitY, 0.0)
>
> earthR :: (Distance, Distance, Distance)
> earthR = (earthPerihelion, 0.0, 0.0)

For completeness we give the Sun's starting conditions.

> sunV :: (Speed, Speed, Speed)
> sunV = (0.0, 0.0, 0.0)
>
> sunR :: (Distance, Distance, Distance)
> sunR = (0.0, 0.0, 0.0)

> initVs :: Array U DIM2 Speed
> initVs = fromListUnboxed (Z :. nBodies :. spaceDim) $ concat xs
>   where
>     nBodies = length xs
>     xs = [ [earthX,   earthY,   earthZ]
>          , [jupiterX, jupiterY, jupiterZ]
>          , [sunX,     sunY,     sunZ]
>          ]
>     (earthX,   earthY,   earthZ)   = earthV
>     (jupiterX, jupiterY, jupiterZ) = jupiterV
>     (sunX,     sunY,     sunZ)     = sunV

> initPsM :: Monad m => m (Array U DIM2 Double)
> initPsM = computeP $ ms2 *^ initVs
>   where
>     (Z :. i :. j) = extent initVs
>     ms2 = extend (Any :. j) masses

> initRs :: Array U DIM2 Distance
> initRs = fromListUnboxed (Z :. nBodies :. spaceDim) $ concat xs
>   where
>     nBodies = length xs
>     xs = [ [earthX,   earthY,   earthZ]
>          , [jupiterX, jupiterY, jupiterZ]
>          , [sunX,     sunY,     sunZ]
>          ]
>     (earthX,   earthY,   earthZ)   = earthR
>     (jupiterX, jupiterY, jupiterZ) = jupiterR
>     (sunX,     sunY,     sunZ)     = sunR

> masses :: Array U DIM1 Mass
> masses = fromListUnboxed (Z :. nBodies) xs
>   where
>     nBodies = length xs
>     xs = [ earthMass
>          , jupiterMass
>          , sunMass
>          ]

> stepN :: forall m . Monad m =>
>          Int -> Masses -> Positions -> Momenta ->
>          m (Positions, Momenta)
> stepN n masses = curry updaterMulti
>   where
>     updaterMulti = foldr (>=>) return updaters
>     updaters :: [(Positions, Momenta) -> m (Positions, Momenta)]
>     updaters = replicate n (uncurry (stepOnceP dt masses))

FIXME: Surely this can be as an instance of some nice recursion pattern.

> stepN' :: Monad m =>
>           Int -> Masses -> Positions -> Momenta ->
>           m [(Positions, Momenta)]
> stepN' n ms rs vs = do
>   rsVs <- stepAux n rs vs
>   return $ (rs, vs) : rsVs
>   where
>     stepAux 0  _  _ = return []
>     stepAux n rs vs = do
>       (newRs, newVs) <- stepOnceP dt ms rs vs
>       rsVs <- stepAux (n-1) newRs newVs
>       return $ (newRs, newVs) : rsVs

Performance
-----------

> nSteps = 36 * 12

> main :: IO ()
> main = do
>   initPs <- initPsM
>   rsVs <- stepN' nSteps masses initRs initPs
>   -- putStrLn $ show rsVs
>   let rs = Prelude.map fst rsVs
>   -- putStrLn $ show $ extent $ rs!!nSteps
>   let vs = Prelude.map snd rsVs
>   putStrLn $ show $ extent $ vs!!nSteps
>   let erx nSteps = (rs!!nSteps)!(Z :. (0 :: Int) :. (0 :: Int))
>       ery nSteps = (rs!!nSteps)!(Z :. (0 :: Int) :. (1 :: Int))
>       erz nSteps = (rs!!nSteps)!(Z :. (0 :: Int) :. (2 :: Int))
>       jrx nSteps = (rs!!nSteps)!(Z :. (1 :: Int) :. (0 :: Int))
>       jry nSteps = (rs!!nSteps)!(Z :. (1 :: Int) :. (1 :: Int))
>       jrz nSteps = (rs!!nSteps)!(Z :. (1 :: Int) :. (2 :: Int))
>       -- srx = (rs!!nSteps)!(Z :. (2 :: Int) :. (0 :: Int))
>       -- sry = (rs!!nSteps)!(Z :. (2 :: Int) :. (1 :: Int))
>       -- srz = (rs!!nSteps)!(Z :. (2 :: Int) :. (2 :: Int))
>   mapM (\n -> putStrLn $ printf "%16.10e %16.10e %16.10e" (erx n) (ery n) (erz n)) [0 .. nSteps - 1]
>   mapM (\n -> putStrLn $ printf "%16.10e %16.10e %16.10e" (jrx n) (jry n) (jrz n)) [0 .. nSteps - 1]
>   -- putStrLn $ printf "%16.10e %16.10e %16.10e" srx sry srz
>   h <- zipWithM (hamiltonianP masses) vs rs
>   putStrLn $ show $ h

> simPlanets = runIdentity $ do
>   initPs <- initPsM
>   rsVs <- stepN' nSteps masses initRs initPs
>   let ps = Prelude.map fst rsVs
>       exs = Prelude.map (!(Z :. (0 :: Int) :. (0 :: Int))) ps
>       eys = Prelude.map (!(Z :. (0 :: Int) :. (1 :: Int))) ps
>       jxs = Prelude.map (!(Z :. (1 :: Int) :. (0 :: Int))) ps
>       jys = Prelude.map (!(Z :. (1 :: Int) :. (1 :: Int))) ps
>       sxs = Prelude.map (!(Z :. (2 :: Int) :. (0 :: Int))) ps
>       sys = Prelude.map (!(Z :. (2 :: Int) :. (1 :: Int))) ps
>   return $ zip3 (zip exs eys) (zip jxs jys) (zip sxs sys)

```{.dia width='600'}
import Symplectic
import SymplecticDia

dia' :: DiagramC

dia' = test tickSize [ (cellColour1, map (\(x, _, _) -> x) simPlanets)
                     , (cellColour2, map (\(_, y, _) -> y) simPlanets)
                     , (cellColour3, map (\(_, _, z) -> z) simPlanets)
                     ]

dia = dia'
```

Bibliography
------------
