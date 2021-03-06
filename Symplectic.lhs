% Planetary Simulation with Excursions in Symplectic Manifolds
% Dominic Steinitz
% 7th June 2013

This article attempts to show that Haskell [@Hudak:2007:HHL:1238844.1238856]
performs reasonably well on numerical problems.

When I started to do this, it seemed straightforward enough: pick a
problem which admitted a numerical solution, find an algorithm and
code it up. I chose the problem of orbital dynamics as I had always
been fascinated by the precession of the perihelion of Mercury (which
is mainly caused by the pull of the other planets in the Solar System)
and because this admits of at least two different methods of numerical
solution both of which I hope will show the power of Haskell in this
area. This led to the selection of a suitable algorithm and I read
that one should prefer a symplectic method such as the Leapfrog which
conserves the energy of a system (a highly desirable requirement when
modelling orbital dynamics). My conscience would not let me write
about such a method without being able to explain it. This led into
the Hamiltonian formulation of classical mechanics, symplectic
manifolds and symplectic (numerical) methods.

The reader interested in the Haskell implementations and performance
comparisons (currently *not* with other programming languages) can
read the [introduction](#Introduction) and skip to the section on
[performance](#Performance). I apologise in advance to experts in
classical mechanics, symplectic geometery and numerical analysis and
can only hope I have not traduced their subjects too much.

Note that we do not make it as far the perihelion of Mercury in this
article but we do simulate the planets in the outer solar system.

Introduction
============

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
\mathbb{H} = \frac{1}{2}ml^2\dot{\theta}^2 + mgl(1 - \cos\theta)
$$

Using the Langrangian ${\mathbb{L}} = T - V$ where $T$ and $V$ are the
kinetic and potential energies respectively, let us set the
generalized momentum

$$
p = \frac{\partial\mathbb{L}}{\partial\dot{\theta}} = ml^2\dot{\theta}
$$

Then we can re-write the Hamiltonian as:

$$
\mathbb{H} = \frac{p^2}{2ml^2} + mgl(1 - \cos\theta)
$$

Applying Hamilton's equations we obtain

$$
\begin{aligned}
\dot{\theta} &=  \frac{\partial\mathbb{H}}{\partial p}      = \frac{p}{ml^2} \\
\dot{p}      &= -\frac{\partial\mathbb{H}}{\partial \theta} = -mgl\sin\theta
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

Haskell for Explicit Euler
--------------------------

First we need some pragmas, exports (required to create the diagrams)
and imports.

> {-# OPTIONS_GHC -Wall                     #-}
> {-# OPTIONS_GHC -fno-warn-name-shadowing  #-}
> {-# OPTIONS_GHC -fno-warn-type-defaults   #-}

> {-# LANGUAGE NoMonomorphismRestriction    #-}
> {-# LANGUAGE FlexibleContexts             #-}
> {-# LANGUAGE ScopedTypeVariables          #-}
> {-# LANGUAGE GeneralizedNewtypeDeriving   #-}
> {-# LANGUAGE TypeOperators                #-}

> module Symplectic (
>     pendulumEE
>   , pendulumSE
>   , jupiterEarth
>   , outerPlanets
>   , main
>   ) where

> import           Data.Array.Repa hiding ((++), zipWith)
> import qualified Data.Array.Repa as Repa

> import           Control.Monad
> import           Control.Monad.Identity
> import           System.Environment
> import           System.Console.GetOpt

> import           Foreign.Storable

> import qualified Data.Yarr as Y
> import           Data.Yarr (loadS, dzip2, dzip3, F, L)
> import           Data.Yarr.Repr.Delayed (UArray)
> import qualified Data.Yarr.Shape as S
> import qualified Data.Yarr.Utils.FixedVector as V
> import           Data.Yarr.Utils.FixedVector (VecList, N3)
> import qualified Data.Yarr.IO.List as YIO
> import qualified Data.Yarr.Walk as W
>
> import qualified Initial as I

Some type synomyms although it is doubtful how useful these are since
the generalized co-ordinates and momenta that one uses with
Hamiltonian methods can have different units depending on how the
physical problem is formulated.

> type Distance = Double
> type Mass     = Double
> type Speed    = Double

The functions to update the position and momentum.

> stepMomentumEE :: Double -> Double -> Double -> Double -> Double
> stepMomentumEE m l p q = p -  h * m * g * l * sin q

> stepPositionEE :: Double -> Double -> Double -> Double -> Double
> stepPositionEE m l p q = q + h * p / (m * l^2)

The explicit Euler method itself. Notice that both update functions
use the previous position and momentum.

> stepOnceEE :: Double -> Double -> Double -> Double -> (Double, Double)
> stepOnceEE m l p q = (newP, newQ)
>   where
>     newP = stepMomentumEE m l p q
>     newQ = stepPositionEE m l p q

The physical data for our problem and also the step length for the numerical method.

> h, m, l, g :: Double
> h = 0.01 -- Seconds
> l = 1.0  -- Metres
> m = 1.0  -- Kilograms
> g = 9.81 -- Metres * Seconds^-2

Let's start our pendulum at the bottom with an angular velocity that
ensures we don't go over the top.

> initTheta, initThetaDot, initP :: Double
> initTheta    = 0.0
> initThetaDot = 1.7
> initP        = m * l^2 * initThetaDot

> runEE :: Double -> Double -> [(Double, Double)]
> runEE initP initTheta = iterate (uncurry (stepOnceEE m l))
>                                 (initP, initTheta)

> pendulumEE :: [(Double, Double)]
> pendulumEE = runEE initP initTheta

The diagram below plots the position of the pendulum (the angle it
makes with the vertical) against momentum, both axes normalised so
that the maximum position and momentum are 1.0. We would expect that
the trajectory would form a closed path that is traversed indefinitely
as the pendulum swings back and forth. Instead we see that trajectory
gradually spirals outward showing that energy is not conserved but
steadily increases over time, an undesirable state of affairs.

```{.dia width='800'}
import Symplectic
import SymplecticDia

diaEE :: DiagramC
diaEE = test tickSize [ (cellColourEE0, take nPlotPoints $ pendulumEE)
                      ]
dia = diaEE
```

Haskell for Symplectic Euler
----------------------------

Instead let us apply the the symplectic Euler method:

$$
\begin{aligned}
p_{n+1} = p_n - hmgl\sin\theta_n \\
\theta_{n+1} = \theta_n + \frac{hp_{n+1}}{2ml^2}
\end{aligned}
$$

The functions to update the position and momentum.

> stepMomentum :: Double -> Double -> Double -> Double -> Double
> stepMomentum m l p q = p -  h * m * g * l * sin q

> stepPosition :: Double -> Double -> Double -> Double -> Double
> stepPosition m l p q = q + h * p / (m * l^2)

The symplectic Euler method itself. Notice that only the update
function for momentum uses both the previous position and momentum;
the update function for position uses the previous position but the
new momentum.

> stepOnce :: Double -> Double -> Double -> Double -> (Double, Double)
> stepOnce m l p q = (newP, newQ)
>   where
>     newP = stepMomentum m l p q
>     newQ = stepPosition m l newP q

> runSE :: Double -> Double -> [(Double, Double)]
> runSE initP initTheta = iterate (uncurry (stepOnce m l))
>                                 (initP, initTheta)

> pendulumSE :: [(Double, Double)]
> pendulumSE = runSE initP initTheta

The diagram below plots the position of the pendulum (the angle it
makes with the vertical) against momentum, both axes normalised so
that the maximum position and momentum are 1.0. We would expect that
the trajectory would form a closed path that is traversed indefinitely
as the pendulum swings back and forth. And indeed this is the case.

```{.dia width='800'}
import Symplectic
import SymplecticDia

dia' :: DiagramC
dia' = test tickSize [ (cellColour0, take nPlotPoints $ pendulumSE)
                     ]

dia = dia'
```

So in this case the energy is conserved so this looks like a good
candidate for simulating orbital dynamics. But why does this work? It
really looks very similar to the explicit Euler method.

Theory
======

Warning: some handwaving as a full and formal exposition of the theory
would take much more space than can be contained in this blog article.

We can think of the evolution of the pendulum as taking place on a
2-dimensional manifold
[manifold](http://en.wikipedia.org/wiki/Manifold "Wikipedia
definition") $\mathbb{S}^1 \times \mathbb{R}$ where $\mathbb{S}^1$ is
the 1-dimensional sphere (a circle) since the pendulum's space
co-ordinate can only take on values $0 \le q < 2\pi$.

We can define a
([symplectic](https://en.wikipedia.org/wiki/Symplectic_manifold
"Wikipedia definition"))
[2-form](https://en.wikipedia.org/wiki/Differential_form "Wikipedia
definition") on this manifold:

$$
\omega = dq \wedge dp
$$

Using this we can now produce a vector field from the Hamiltonian:
$\mathbb{H} : \mathbb{S}^1 \times \mathbb{R} \longrightarrow \mathbb{R}$

In order to this and without proof let us record the following fact.

**Theorem**

Let $(M, \omega)$ be a symplectic manifold. Then there exists a
bundle isomorphism $\tilde{\omega} : TM \longrightarrow T^*M$ defined
by $\tilde{\omega}(X_p)(Y_p) = \omega_p(X_p, Y_p)$.

$\blacksquare$

This is analagous to the isomorphism one can derive in a (semi)
Riemannian manifold with the metric in some sense playing the role of
the 2-form (see [@o1983semi] for example).

We assume the Hamiltonian to be a smooth function. We can form the
1-form $dH$ and we can define the Hamiltonian vector field $X_H =
\tilde{\omega}^{-1}(dH)$.

We have

$$
d\mathbb{H} = \frac{\partial{\mathbb{H}}}{\partial q}dq +
           \frac{\partial{\mathbb{H}}}{\partial p}dp
$$

Thus the corresponding vector field is given by

$$
X_\mathbb{H} =  \frac{\partial{\mathbb{H}}}{\partial q}\frac{\partial}{\partial q} -
             \frac{\partial{\mathbb{H}}}{\partial p}\frac{\partial}{\partial p}
$$

The flow of this vector field is the solution to

$$
\begin{aligned}
\dot{q} &=  \frac{\partial \mathbb{H}}{\partial p} \\
\dot{p} &= -\frac{\partial \mathbb{H}}{\partial q} \\
\end{aligned}
$$

In other words by using the symplectic 2-form and the Hamiltonian we
have regained Hamilton's equations.

**Theorem**

*$\mathbb{H}$ is constant on flows of $X_\mathbb{H}$.*

*Proof*

$$
X_{\mathbb{H}}{\mathbb{H}} = \omega(X_{\mathbb{H}}, X_{\mathbb{H}}) = 0
$$

since $\omega$ is alternating.

$\blacksquare$

When the Hamiltonian function represents the energy of the system
being studied then this says that energy remains constant on
flows. That is, as the system evolves according to the flow $\phi_t$
given by the vector field $X_{\mathbb{H}}$ then $\mathbb{H}(q_t, p_t) =
\mathbb{H}(\phi_t(q_0, p_0)) = \mathbb{H}(q_0, p_0)$.

Thus it makes sense to look for numeric methods which maintain this
invariant, that is methods which preserve the symplectic 2-form.

**Definition**

A diffeomorphsim between two symplectic manifolds $f : (M, \mu)
\longrightarrow (M, \nu)$ is a symplectomorphism if

$$
f^*\nu = \mu
$$

$\blacksquare$

In co-ordinates, we have

$$
\begin{aligned}
f^*(dq \wedge dp) &= (\frac{\partial f_q}{\partial q} dq + \frac{\partial f_q}{\partial p} dp)
                     \wedge
                     (\frac{\partial f_p}{\partial q} dq + \frac{\partial f_p}{\partial p} dp) \\
                  &= (\frac{\partial f_q}{\partial q}\frac{\partial f_p}{\partial p} -
                      \frac{\partial f_p}{\partial q}\frac{\partial f_q}{\partial p})dq \wedge dp \\
                  &= dq \wedge dp
\end{aligned}
$$

Or in matrix form

$$
{\begin{bmatrix}
\frac{\partial f_q}{\partial q} & \frac{\partial f_q}{\partial p} \\
\frac{\partial f_p}{\partial q} & \frac{\partial f_p}{\partial p}
\end{bmatrix}}^\top
\,
\begin{bmatrix}
 0 & 1 \\
-1 & 0
\end{bmatrix}
\,
\begin{bmatrix}
\frac{\partial f_q}{\partial q} & \frac{\partial f_q}{\partial p} \\
\frac{\partial f_p}{\partial q} & \frac{\partial f_p}{\partial p}
\end{bmatrix}
=
\begin{bmatrix}
 0 & 1 \\
-1 & 0
\end{bmatrix}
$$

We now check that the symplectic Euler method satisfies this.

Symplectic Euler
----------------

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

Now substitute in the functions for the Euler symplectic method and we
obtain

$$
\begin{bmatrix}
1 + hH_{qp}  & 0 \\
-hG_{pp} & 1
\end{bmatrix}
\,
\frac{\partial \big(p_{n+1},q_{n+1}\big)}{\partial \big(p_{n},q_{n}\big)}
=
\begin{bmatrix}
1 & -hH_{qq} \\
0 & 1 + hH_{qp}
\end{bmatrix}
$$

The reader can then check by a straightforward but tedious calculation
that

$$
\frac{\partial \big(p_{n+1},q_{n+1}\big)}{\partial \big(p_{n},q_{n}\big)}^\top J \frac{\partial \big(p_{n+1},q_{n+1}\big)}{\partial \big(p_{n},q_{n}\big)} = J
$$

where

$$
J
=
\begin{bmatrix}
 0 & 1 \\
-1 & 0
\end{bmatrix}
$$

Thus the symplectic Euler method really is symplectic.

\blacksquare

On the other hand for the explicit Euler for this particular example
we have

$$
\frac{\partial \big(p_{n+1},q_{n+1}\big)}{\partial \big(p_{n},q_{n}\big)}
 =
\begin{bmatrix}
 1 & h / ml^2 \\
-hmglcos\theta_n & 1
\end{bmatrix}
$$

And a simple calculation shows that

$$
\frac{\partial \big(p_{n+1},q_{n+1}\big)}{\partial \big(p_{n},q_{n}\big)}^\top J \frac{\partial \big(p_{n+1},q_{n+1}\big)}{\partial \big(p_{n},q_{n}\big)}
=
\begin{bmatrix}
 0 & 1 + h^2\cos\theta\\
-1 - h^2\cos\theta & 0
\end{bmatrix}
$$

Thus the explicit Euler method is not symplectic i.e., does not
preserve areas. Thus the path traversed is not an integral curve of
the Hamiltonian vector field. We can see this in the diagram: the path
spirals outwards. More details and examples can be found in
[@IAUS:152:407Y; @Cross:2005:SIOC].

Planetary Motion
================

Normally we would express the gravitational constant in SI units but
to be consistent with [@hairer2010geometric] we use units in which
distances are expressed in astronomical units, masses are measured
relative to the sun and time is measured in earth days.

$$
{\mathbb H} = \frac{1}{2}\sum_{i=0}^n \frac{p_i^\top p_i}{m_i} - \frac{G}{2}\sum_{i=0}^n\sum_{j \neq i} \frac{m_i m_j}{\|q_i - q_j\|}
$$

Applying Hamilton's equations we obtain

$$
\begin{aligned}
\dot{q_k^a} &=  \frac{\partial\mathbb{H}}{\partial p_k^a} = \frac{p_k^a}{m_k} \\
\dot{p_k^a} &= -\frac{\partial\mathbb{H}}{\partial q_k^a} = G\sum_{j \neq k}m_k m_i \frac{q_k^a - q_j^a}{\|q_k - q_j\|^3}
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

Repa Implementation
-------------------

We use [repa](http://hackage.haskell.org/package/repa "Hackage")
represent our positions and momenta as 2-dimensional arrays, each
planet is given a 3-dimensional position vector and a 3-dimensional
momentum vector.

> newtype PositionP a = QP { positionP :: Array a DIM2 Double }
> newtype MomentaP  a = PP { momentaP :: Array a DIM2 Double }
> newtype MassP     a = MP { massP :: Array a DIM1 Double }

> stepPositionP :: forall a b c m .
>                  ( Monad m
>                  , Source a Double
>                  , Source b Double
>                  , Source c Double
>                  ) =>
>                  Double ->
>                  PositionP a ->
>                  MassP b ->
>                  MomentaP c ->
>                  m (PositionP U)
> stepPositionP h qs ms ps = do
>   do newQs <- computeP $ (positionP qs) +^ ((momentaP ps) *^ h2 /^ ms2)
>      return $ QP newQs
>     where
>       (Z :. i :. j) = extent $ momentaP ps
>
>       h2  = extend (Any :. i :. j) $ fromListUnboxed Z [h]
>       ms2 = extend (Any :. j) $ massP ms

Each planet produces forces on every other planet so we work with
3-dimsenional arrays and explicitly set the force of a planet on
itself to zero. Once the forces on each planet have been calculated,
we sum them to produce a resultant force which we then use to step the
momentum forward.

> stepMomentumP :: forall a b c m .
>                  ( Monad m
>                  , Source a Double
>                  , Source b Double
>                  , Source c Double
>                  ) =>
>                  Double ->
>                  Double ->
>                  PositionP a ->
>                  MassP b ->
>                  MomentaP c ->
>                  m (MomentaP U)
> stepMomentumP gConst h qs ms ps =
>   do fs <- sumP $ transpose $ zeroDiags fss
>      newPs <- computeP $ (momentaP ps) +^ (fs *^ dt2)
>      return $ PP newPs
>   where
>     is = repDim2to3Outer $ prodPairsMasses $ massP ms
>     qDiffs = pointDiffs $ positionP qs
>     preDs = Repa.map (^3) $
>             Repa.map sqrt $
>             sumS $
>             Repa.map (^2) $
>             qDiffs
>     ds    = repDim2to3Outer preDs
>     preFs = Repa.map (* (negate gConst)) $
>             qDiffs /^ ds
>     fss = is *^ preFs
>     Z :.i :. _j :. k = extent fss
>     dt2              = extend (Any :. i :. k) $ fromListUnboxed Z [h]
>
>     repDim2to3Outer a = extend (Any :. I.spaceDim) a
>
>     zeroDiags x = traverse x id f
>       where
>         f _ (Z :. i :. j :. k) | i == j    = 0.0
>                                | otherwise = x!(Z :. i :. j :. k)

> stepOnceP :: ( Monad m
>              , Source a Double
>              , Source b Double
>              , Source c Double
>              ) =>
>              Double ->
>              Double ->
>              MassP b ->
>              PositionP a ->
>              MomentaP c ->
>              m (PositionP U, MomentaP U)
> stepOnceP gConst h ms qs ps = do
>   newPs <- stepMomentumP gConst h qs ms ps
>   newQs <- stepPositionP h qs ms newPs
>   return (newQs, newPs)


> kineticEnergyP :: MassP U -> MomentaP U -> IO (Array D DIM0 Double)
> kineticEnergyP ms ps = do
>   preKes <- sumP $ (momentaP ps) *^ (momentaP ps)
>   ke     <- sumP $ preKes /^ (massP ms)
>   return $ Repa.map (* 0.5) ke
>
> potentialEnergyP :: Double ->
>                     MassP U ->
>                     PositionP U ->
>                     IO (Array U DIM1 Double)
> potentialEnergyP gConst ms qs = do
>   ds2 <- sumP $ Repa.map (^2) $ pointDiffs $ positionP qs
>   let ds   = Repa.map sqrt ds2
>       is   = prodPairsMasses $ massP ms
>       pess = zeroDiags $ Repa.map (* (0.5 * negate gConst)) $ is /^ ds
>   pes <- sumP pess
>   return pes
>
>   where
>
>     zeroDiags x = traverse x id f
>       where
>         f _ (Z :. i :. j) | i == j    = 0.0
>                           | otherwise = x!(Z :. i :. j)

> hamiltonianP :: Double ->
>                 MassP U ->
>                 PositionP U ->
>                 MomentaP U ->
>                 IO Double
> hamiltonianP gConst ms qs ps = do
>   ke <- kineticEnergyP ms ps
>   pes <- potentialEnergyP gConst ms qs
>   pe  <- sumP pes
>   te :: Array U DIM0 Double <- computeP $ ke +^ pe
>   return $ head $ toList te

> prodPairsMasses :: Source a Double =>
>                    Array a DIM1 Double ->
>                    Array D DIM2 Double
> prodPairsMasses ms = ns *^ (transpose ns)
>
>   where
>     (Z :. i) = extent ms
>     ns = extend (Any :. i :. All) ms

>
> pointDiffs :: Source a Double =>
>               Array a DIM2 Double ->
>               Array D DIM3 Double
> pointDiffs qs = qss -^ (transposeOuter qss)
>   where
>
>     qss = replicateRows qs
>
>     transposeOuter qs = backpermute (f e) f qs
>       where
>         e = extent qs
>         f (Z :. i :. i' :. j) = Z :. i' :. i :. j
>
>     replicateRows :: Source a Double =>
>                      Array a DIM2 Double ->
>                      Array D DIM3 Double
>     replicateRows a = extend (Any :. i :. All) a
>       where (Z :. i :. _j) = extent a

Using the single step udate, we can step as many times as we wish.

> stepN :: forall m . Monad m =>
>          Int ->
>          Double ->
>          Double ->
>          MassP U ->
>          PositionP U ->
>          MomentaP U ->
>          m (PositionP U, MomentaP U)
> stepN n gConst dt masses = curry updaterMulti
>   where
>     updaterMulti = foldr (>=>) return updaters
>     updaters = replicate n (uncurry (stepOnceP gConst dt masses))

Sometimes we need all the intermediate steps e.g. for plotting.

> stepNs :: Monad m =>
>           Int ->
>           Double ->
>           Double ->
>           MassP U ->
>           PositionP U ->
>           MomentaP U ->
>           m [(PositionP U, MomentaP U)]
> stepNs n gConst dt ms rs vs = do
>   rsVs <- stepAux n rs vs
>   return $ (rs, vs) : rsVs
>   where
>     stepAux 0  _  _ = return []
>     stepAux n rs vs = do
>       (newRs, newVs) <- stepOnceP gConst dt ms rs vs
>       rsVs <- stepAux (n-1) newRs newVs
>       return $ (newRs, newVs) : rsVs

Yarr Implementation
-------------------

We use [yarr](http://hackage.haskell.org/package/yarr "Hackage")
represent our positions and momenta as 1-dimensional arrays, each
planet is given a 3-dimensional position vector and a 3-dimensional
momentum vector.

> vZero :: VecList N3 Double
> vZero = V.replicate 0

> type ArrayY = UArray F L S.Dim1

> newtype PositionY  = QY { positionY :: VecList N3 Double }
>   deriving (Show, Storable)
> newtype MomentumY = PY { momentumY :: VecList N3 Double }
>   deriving (Show, Storable)
> type MomentaY     = ArrayY MomentumY
> type PositionsY   = ArrayY PositionY
> type MassesY      = ArrayY Mass
> type ForceY       = VecList N3 Double
> type ForcesY      = ArrayY ForceY

> stepPositionY :: Double -> PositionsY -> MassesY -> MomentaY -> IO ()
> stepPositionY h qs ms ps = do
>   loadS S.fill (dzip3 upd qs ms ps) qs
>   where
>     upd :: PositionY -> Mass -> MomentumY -> PositionY
>     upd q m p = QY $ V.zipWith (+)
>                 (positionY q)
>                 (V.map (* (h / m)) (momentumY p))

Note the requirement to fill the forces array with zeros before using it.

> stepMomentumY :: Double ->
>                  Double ->
>                  PositionsY ->
>                  MassesY ->
>                  MomentaY ->
>                  IO ()
> stepMomentumY gConst h qs ms ps = do
>   let nBodies = Y.extent ms
>   fs :: ForcesY <- Y.new nBodies
>   S.fill (\_ -> return vZero) (Y.write fs) 0 nBodies
>   let forceBetween i pos1 mass1 j
>         | i == j = return vZero
>         | otherwise = do
>           pos2 <- qs `Y.index` j
>           mass2 <- ms `Y.index` j
>           let deltas = V.zipWith (-) (positionY pos1) (positionY pos2)
>               dist2  = V.sum $ V.map (^ 2) deltas
>               a = 1.0 / dist2
>               b = (negate gConst) * mass1 * mass2 * a * (sqrt a)
>           return $ V.map (* b) deltas
>       forceAdd :: Int -> Int -> ForceY -> IO ()
>       forceAdd i _ f = do
>         f0 <- fs `Y.index` i
>         Y.write fs i (V.zipWith (+) f0 f)
>       force i pos = do
>         mass <- ms `Y.index` i
>         S.fill (forceBetween i pos mass) (forceAdd i) 0 nBodies
>       upd momentum force =
>         PY $ V.zipWith (+)
>         (momentumY momentum)
>         (V.map (\f -> f * h) force)
>   S.fill (Y.index qs) force 0 nBodies
>   loadS S.fill (dzip2 upd ps fs) ps

> stepOnceY :: Double ->
>              Double ->
>              MassesY ->
>              PositionsY ->
>              MomentaY ->
>              IO ()
> stepOnceY gConst h ms qs ps = do
>   stepMomentumY gConst h qs ms ps
>   stepPositionY h qs ms ps

> potentialEnergyY :: Double ->
>                     MassesY ->
>                     PositionsY ->
>                     IO (ArrayY Double)
> potentialEnergyY gConst ms qs = do
>   let nBodies = Y.extent ms
>   pes :: ArrayY Double <- Y.new nBodies
>   S.fill (\_ -> return 0.0) (Y.write pes) 0 nBodies
>   let peOnePairParticles :: Int ->
>                             Int ->
>                             IO Double
>       peOnePairParticles i j
>         | i == j = return 0.0
>         | otherwise = do
>           q1 <- qs `Y.index` i
>           m1 <- ms `Y.index` i
>           q2 <- qs `Y.index` j
>           m2 <- ms `Y.index` j
>           let qDiffs = V.zipWith (-) (positionY q1) (positionY q2)
>               dist2  = V.sum $ V.map (^2) qDiffs
>               a      = 1.0 / dist2
>               b      = 0.5 * (negate gConst) * m1 * m2 * (sqrt a)
>           return b
>       peAdd i _ pe = do
>         peDelta <- pes `Y.index` i
>         Y.write pes i (pe + peDelta)
>       peFn i _ = do
>         S.fill (peOnePairParticles i) (peAdd i) 0 nBodies
>   S.fill (Y.index qs) peFn 0 nBodies
>   return pes

> kineticEnergyY :: MassesY -> MomentaY-> IO Double
> kineticEnergyY ms ps = do
>   let nakedPs = Y.delay $ Y.dmap momentumY ps
>   let preKes = Y.dmap V.sum $ dzip2 (V.zipWith (*)) nakedPs nakedPs
>       kes     = dzip2 (/) preKes (Y.delay ms)
>   ke <- W.walk (W.reduceL S.foldl (+)) (return 0) kes
>   return $ 0.5 * ke

> hamiltonianY :: Double -> MassesY -> PositionsY -> MomentaY-> IO Double
> hamiltonianY gConst ms qs ps = do
>   ke  <- kineticEnergyY ms ps
>   pes <- potentialEnergyY gConst ms qs
>   pe  <- W.walk (W.reduceL S.foldl (+)) (return 0) pes
>   return $ pe + ke

The Outer Solar System
======================

We convert the starting positions and momenta for the planets into
repa and yarr friendly representations.

> mosssP :: MassP U
> mosssP = MP $ fromListUnboxed (Z :. n) I.massesOuter
>   where
>     n = length I.massesOuter
>
> mosssY :: IO MassesY
> mosssY = YIO.fromList n I.massesOuter
>   where
>     n = length I.massesOuter
>
> qosss :: PositionP U
> qosss = QP $ fromListUnboxed (Z :. n :. I.spaceDim) xs
>   where
>     xs = concat I.initQsOuter
>     n  = length xs `div` I.spaceDim
>
> qosssY :: IO PositionsY
> qosssY = YIO.fromList nBodies $ Prelude.map f [0 .. nBodies - 1]
>   where
>     nBodies = length I.initQsOuter
>     f :: Int -> PositionY
>     f i = QY $
>           V.vl_3 ((I.initQsOuter!!i)!!0)
>                  ((I.initQsOuter!!i)!!1)
>                  ((I.initQsOuter!!i)!!2)
>
> posss :: MomentaP U
> posss = PP $ fromListUnboxed (Z :. n :. I.spaceDim) xs
>   where
>     xs = concat I.initPsOuter
>     n  = length xs `div` I.spaceDim
>
> posssY :: IO MomentaY
> posssY = YIO.fromList nBodies $ Prelude.map f [0 .. nBodies - 1]
>   where
>     nBodies = length I.initPsOuter
>     f :: Int -> MomentumY
>     f i = PY $
>           V.vl_3 ((I.initPsOuter!!i)!!0)
>                  ((I.initPsOuter!!i)!!1)
>                  ((I.initPsOuter!!i)!!2)

Rather arbitrarily we run the outer planets for 2000 steps with a step
length of 100 days.

> outerPlanets :: [[(Double, Double)]]
> outerPlanets = runIdentity $ do
>   rsVs <- stepNs 2000 I.gConstAu 100 mosssP qosss posss
>   let qs = Prelude.map fst rsVs
>       xxs = Prelude.map
>             (\i -> Prelude.map ((!(Z :. (i :: Int) :. (0 :: Int))) .
>                                 positionP) qs)
>             [5,0,1,2,3,4]
>       xys = Prelude.map
>             (\i -> Prelude.map ((!(Z :. (i :: Int) :. (1 :: Int))) .
>                                 positionP) qs)
>             [5,0,1,2,3,4]
>   return $ zipWith zip xxs xys

Plotting the results, we can see that we have a simulation which, as
we expect, conserves energy.

```{.dia width='600'}
import Symplectic
import SymplecticDia

dia' :: DiagramC

dia' = test tickSize [ (cellColour0, zipWith (-) (outerPlanets!!0) (outerPlanets!!0))
                     , (cellColour1, zipWith (-) (outerPlanets!!1) (outerPlanets!!0))
                     , (cellColour2, zipWith (-) (outerPlanets!!2) (outerPlanets!!0))
                     , (cellColour1, zipWith (-) (outerPlanets!!3) (outerPlanets!!0))
                     , (cellColour2, zipWith (-) (outerPlanets!!4) (outerPlanets!!0))
                     , (cellColour3, zipWith (-) (outerPlanets!!5) (outerPlanets!!0))
                     ]

dia = dia'
```

Performance
===========

Let's see how repa and yarr perform against each other.

> data YarrOrRepa = Repa | Yarr
>   deriving Show
>
> data Options = Options  { optYarr :: YarrOrRepa
>                         }
>
> startOptions :: Options
> startOptions = Options  { optYarr = Repa
>                         }
>
> options :: [OptDescr (Options -> IO Options)]
> options = [
>   Option ['Y'] ["yarr"] (NoArg (\opt -> return opt { optYarr = Yarr }))
>          "Use yarr"
>   ]

> main :: IO ()
> main = do
>   args <- getArgs
>   let (actions, _nonOpts, _msgs) = getOpt RequireOrder options args
>   opts <- foldl (>>=) (return startOptions) actions
>   case optYarr opts of
>     Repa -> do
>       hPre <- hamiltonianP I.gConstAu mosssP qosss posss
>       putStrLn $ show hPre
>       (qsPost, psPost) <- stepN I.nStepsOuter I.gConstAu I.stepOuter
>                           mosssP qosss posss
>       hPost <- hamiltonianP I.gConstAu mosssP qsPost psPost
>       putStrLn $ show hPost
>     Yarr -> do
>       ms :: MassesY <- mosssY
>       ps <- posssY
>       qs <- qosssY
>       hPre <- hamiltonianY I.gConstAu ms qs ps
>       putStrLn $ show hPre
>       S.fill (\_ -> return ())
>              (\_ _ -> stepOnceY I.gConstAu I.stepOuter ms qs ps)
>              (0 :: Int) I.nStepsOuter
>       hPost <- hamiltonianY I.gConstAu ms qs ps
>       putStrLn $ show hPost


With 200,000 steps we get the following for repa.

    $ time ./Symplectic
    -3.215453183208164e-8
    -3.139737384661333e-8

    real	0m18.400s
    user	0m18.245s
    sys	0m0.154s

And a significant speed up with yarr.

    $ time ./Symplectic -Y
    -3.215453183208164e-8
    -3.13973738466838e-8

    real	0m0.553s
    user	0m0.539s
    sys	0m0.013s

With 2,000,000 steps (about 548,000 years) we get the following with
yarr.

    $ time ./Symplectic -Y
    -3.215453183208164e-8
    -3.2144315777817145e-8

    real	0m5.477s
    user	0m5.369s
    sys	0m0.107s

It would be interesting to compare this against a C implementation.

Appendices
==========

Appendix A: Jupiter, Earth and Sun
----------------------------------

We need some initial conditions to start our simulation. Instead of
taking real data, let's make up something which is realistic but
within our control.

  [jupiter]: http://en.wikipedia.org/wiki/Jupiter

Following [@Fitz:Newtonian:Dynamics] we can write Kepler's laws as

$$
\begin{aligned}
r &= \frac{a(1 - e^2)}{1 - e\cos\theta} \\
r^2\dot{\theta} &= \sqrt{(1 - e^2)}na^2 \\
GM_{\rm Sun} &= n^2a^3
\end{aligned}
$$

where $T$ is the period of the orbit, $n = 2\pi / T$ is the mean
angular orbital velocity, $a$ is the major radius of the elliptical
orbit (Kepler's first law: "The orbit of every planet is an ellipse
with the Sun at one of the two foci"), $G$ is the gravitational
constant and $M_{\rm Sun}$ is the mass of the sun and $e$ is the
eccentricity of a planet's orbit.

We can calculate the mean angular orbital velocity of Jupiter:

> nJupiter :: Double
> nJupiter = sqrt $ I.gConst * I.sunMass / I.jupiterMajRad^3

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

In Haskell, we get the following initial conditions:

> jupiterThetaDotP :: Double
> jupiterThetaDotP =
>   nJupiter *
>   I.jupiterMajRad^2 *
>   sqrt (1 - I.jupiterEccentrity^2) / I.jupiterPerihelion^2
>
> jupiterDeltaThetaP :: Double
> jupiterDeltaThetaP = jupiterThetaDotP * I.stepTwoPlanets / 2
>
> jupiterVPeri :: Speed
> jupiterVPeri = jupiterThetaDotP * I.jupiterPerihelion
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
> jupiterR = (negate I.jupiterPerihelion, 0.0, 0.0)

We can do the same for Earth but we assume the earth is at its
perihelion on the opposite side of the Sun to Jupiter.

> nEarth :: Double
> nEarth = sqrt $ I.gConst * I.sunMass / I.earthMajRad^3
>
> earthThetaDotP :: Double
> earthThetaDotP = nEarth *
>                  I.earthMajRad^2 *
>                  sqrt (1 - I.earthEccentrity^2) / I.earthPerihelion^2
>
> earthDeltaThetaP :: Double
> earthDeltaThetaP = earthThetaDotP * I.stepTwoPlanets / 2
>
> earthVPeri :: Speed
> earthVPeri = earthThetaDotP * I.earthPerihelion
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
> earthR = (I.earthPerihelion, 0.0, 0.0)

For completeness we give the Sun's starting conditions.

> sunV :: (Speed, Speed, Speed)
> sunV = (0.0, 0.0, 0.0)
>
> sunR :: (Distance, Distance, Distance)
> sunR = (0.0, 0.0, 0.0)

> initVs :: Array U DIM2 Speed
> initVs = fromListUnboxed (Z :. nBodies :. I.spaceDim) $ concat xs
>   where
>     nBodies = length xs
>     xs = [ [earthX,   earthY,   earthZ]
>          , [jupiterX, jupiterY, jupiterZ]
>          , [sunX,     sunY,     sunZ]
>          ]
>     (earthX,   earthY,   earthZ)   = earthV
>     (jupiterX, jupiterY, jupiterZ) = jupiterV
>     (sunX,     sunY,     sunZ)     = sunV

> initPs :: MomentaP U
> initPs = PP $ runIdentity $ computeP $ ms2 *^ initVs
>   where
>     (Z :. _i :. j) = extent initVs
>     ms2 = extend (Any :. j) (massP masses)

> initQs :: PositionP U
> initQs = QP $ fromListUnboxed (Z :. nBodies :. I.spaceDim) $ concat xs
>   where
>     nBodies = length xs
>     xs = [ [earthX,   earthY,   earthZ]
>          , [jupiterX, jupiterY, jupiterZ]
>          , [sunX,     sunY,     sunZ]
>          ]
>     (earthX,   earthY,   earthZ)   = earthR
>     (jupiterX, jupiterY, jupiterZ) = jupiterR
>     (sunX,     sunY,     sunZ)     = sunR

> masses :: MassP U
> masses = MP $ fromListUnboxed (Z :. nBodies) I.massesTwoPlanets
>   where
>     nBodies = length I.massesTwoPlanets

> jupiterEarth :: [((Double, Double),
>                   (Double, Double),
>                   (Double, Double))]
> jupiterEarth = runIdentity $ do
>   rsVs <- stepNs I.nStepsTwoPlanets I.gConst I.stepTwoPlanets
>                  masses initQs initPs
>   let qs = Prelude.map fst rsVs
>       exs = Prelude.map ((!(Z :. (0 :: Int) :. (0 :: Int))) .
>                          positionP) qs
>       eys = Prelude.map ((!(Z :. (0 :: Int) :. (1 :: Int))) .
>                          positionP) qs
>       jxs = Prelude.map ((!(Z :. (1 :: Int) :. (0 :: Int))) .
>                          positionP) qs
>       jys = Prelude.map ((!(Z :. (1 :: Int) :. (1 :: Int))) .
>                          positionP) qs
>       sxs = Prelude.map ((!(Z :. (2 :: Int) :. (0 :: Int))) .
>                          positionP) qs
>       sys = Prelude.map ((!(Z :. (2 :: Int) :. (1 :: Int))) .
>                          positionP) qs
>   return $ zip3 (zip exs eys) (zip jxs jys) (zip sxs sys)

Plotting the results we can see a reasonable picture for Jupiter's and
Earth's orbits.

```{.dia width='600'}
import Symplectic
import SymplecticDia

dia' :: DiagramC

dia' = test tickSize [ (cellColour0, map (\(x, _, _) -> x) jupiterEarth)
                     , (cellColour1, map (\(_, y, _) -> y) jupiterEarth)
                     , (cellColour2, map (\(_, _, z) -> z) jupiterEarth)
                     ]

dia = dia'
```

Appendix B: The Canonical Symplectic Form for the Cotangent Bundle
------------------------------------------------------

There are plenty of symplectic manifolds besides $\mathbb{R}^{2n}$. The
cotangent bundle has a canonical symplectic 2-form and hence is a
symplectic manifold.

*Proof*

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

\blacksquare

Bibliography
------------
