> {-# LANGUAGE NoMonomorphismRestriction #-}
> {-# LANGUAGE FlexibleContexts          #-}
> {-# LANGUAGE TypeOperators             #-}

> {-# OPTIONS_GHC -Wall                    #-}
> {-# OPTIONS_GHC -fno-warn-name-shadowing #-}
> {-# OPTIONS_GHC -fno-warn-type-defaults  #-}

> import Data.Array.Repa as Repa hiding ((++))

> import Data.Random ()
> import Data.Random.Distribution.Uniform
> import Data.RVar

> import System.Random
> import Control.Monad
> import Control.Monad.State

We wish to model planetary motion using the
[leapfrog][leapfrog] method. This is preferred to other numerical methods as it maintains the total energy of the system.

  [leapfrog]: http://en.wikipedia.org/wiki/Leapfrog_integration

In essence, we have to update the position and velocity of each planet
half a time step out of phase (hence the name leapfrog) as shown
below.

$$
\begin{align*}
x_i                 &= x_{i-1} + v_{i - \frac{1}{2}}\Delta t \\
a_i                 &= F(x_i) \\
v_{i + \frac{1}{2}} &= v_{i - \frac{1}{2}} + a_i\Delta t
\end{align*}
$$

We work in a 3 dimensional Euclidean space.

> spaceDim :: Int
> spaceDim = 3

It's easy to confuse units so we define some type synomyms to help. Of
course having these checked by machine would be better.

FIXME: Can we use dimensional
(http://hackage.haskell.org/package/dimensional-0.10.2) to enforce
units? My experience was not good but maybe that was just me. Perhaps
Chris could look at this, we can contact the author and scour the
Internet.

> type Distance = Double
> type Mass     = Double
> type Force    = Double
> type Speed    = Double
> type Time     = Double

Some physical constants for our system.

> gConst, k, r, mass :: Double
> gConst = 6.67384e-11             -- gravitational constant
> k = 24*60*60                -- seconds in a day
> r = 1e1                     -- radius of sphere particles contained in
> mass = 1e24                 -- mass

> nBodies :: Int
> nBodies = 3                -- number of bodies

And some constants for our finite difference method.

> eps, days, t, timestepDays,dt :: Double
> eps = 0.1*r                 -- softening constant
> days = 36500*100            -- total time in days
> t = days*k                  -- total time
> timestepDays = 10           -- timestep in days
> dt = timestepDays*k         -- timestep
>
> nIters :: Int
> nIters = floor (t/dt)       -- number of iterations

We represent our system of bodies of identical mass as a
1-dimensional (unboxed) array.

> m :: Array D DIM1 Double
> m = Repa.map (mass*) $
>     fromListUnboxed (Z :. nBodies) $
>     take nBodies $ repeat 1.0
>
> rands :: Int -> Double -> Double -> [Double]
> rands n a b =
>   fst $ runState (replicateM n (sampleRVar (uniform a b))) (mkStdGen seed)
>     where
>       seed = 0

We assign our bodies random positions inside a box using a
2-dimensional (unboxed) array.

> randomPoss :: Int -> Double -> Array U DIM2 Double
> randomPoss n r = fromListUnboxed (Z :. n :. spaceDim) $
>                  Prelude.map (fromIntegral . round) $
>                  rands (spaceDim * n) 0 r
>
> startPoss :: Array U DIM2 Double
> startPoss = randomPoss nBodies r

We wish to calculate a two dimensional array of forces where each
force is itself a one dimensional vector of three elements (in the
case of the three dimensional Euclidean space in which we are
operating).

$$
F_i = G\sum_{j\neq i} -m_i m_j\frac{\vec{r}_i - \vec{r}_j}{{\|{\vec{r}_i - \vec{r}_j}\|}^3}
$$

We wish to calculate the products of all pairs of masses. To do this
we first take our 1-dimensional array of masses and create a
2-dimensional array by duplicating the 1-dimensional array $n$ times
where $n$ is the size of the 1-dimensional array. We thus create a
2-dimensional array containing $n \times n$ elements.

> repDim1To2Inner :: Source a Double =>
>                    Array a DIM1 Double ->
>                    Array D DIM2 Double
> repDim1To2Inner a = extend (Any :. i :. All) a
>   where (Z :. i) = extent a

We can now multiply this 2-dimensional array with its transpose
pointwise to get all the products of pairs. Note that the resulting
array is symmetric so we have performed roughly twice as many
multiplications as strictly necessary.

> prodPairsMasses :: Source a Double =>
>                    Array a DIM1 Double ->
>                    Array D DIM2 Double
> prodPairsMasses ms = Repa.zipWith (*) ns (transpose ns)
>   where
>     ns = repDim1To2Inner ms

Next we wish to calculate the distances between each point and every
other point. In a similar way to calculating the pairs of products of
the masses we first replicate a 1-dimensional array of positions in
3-dimensional Euclidean space which we represent as a 2-dimensional
array. We thus create a 3-dimensional array constaining $n \times n
\times 3$ elements.

> replicateRows :: Source a Double =>
>                  Array a DIM2 Double ->
>                  Array D DIM3 Double
> replicateRows a = extend (Any :. i :. All) a
>   where (Z :. i :. _j) = extent a

In this case we cannot use tranpose directly as this will tranpose the
two outermost dimensions and we need to transpose the two innermost
dimensions.

> transposeInner :: Source a Double =>
>               Array a DIM3 Double ->
>               Array D DIM3 Double
> transposeInner ps = backpermute (f e) f ps
>   where
>     e = extent ps
>     f (Z :. i :. i' :. j) = Z :. i' :. i :. j

Now can calculate the point differences between all points and all
other points. Note again that this array is symmetric so we have
performed roughly twice as many multiplications as necessary.

> pointDiffs :: Source a Double =>
>               Array a DIM2 Double ->
>               Array D DIM3 Double
> pointDiffs ps = Repa.zipWith (-) qs (transposeInner qs)
>   where qs = replicateRows ps

We have an 3-dimensional array of distances between each point which
is in effect a 2-dimensional array of the $x$, $y$, and $z$
co-ordinates. We need to mulitply these distances by the cross
products of the masses. But the cross products of the masses form a
2-dimensional array. So we replicate the cross products for each of
the $x$, $y$ and $z$ co-ordinates to form a 3-dimensional array.

> repDim2to3Outer :: Source a Double =>
>              Array a DIM2 Double -> Array D DIM3 Double
> repDim2to3Outer a = extend (Any :. spaceDim) a
>

Now we can calculate the forces in repa using the above equation.

> forces :: (Source a Double, Source b Double) =>
>           Array a DIM2 Double ->
>           Array b DIM1 Double ->
>           Array D DIM3 Double
> forces qs ms = Repa.zipWith (*) fs is
>   where ds = repDim2to3Outer $
>              Repa.map sqrt $
>              Repa.map (+ (eps^2)) $
>              foldS (+) 0 $
>              Repa.map (^2) $
>              pointDiffs qs
>         is = repDim2to3Outer $ prodPairsMasses ms
>         fs = Repa.map (* (negate gConst)) $
>              Repa.zipWith (/) (pointDiffs qs) ds

Next we turn our attention to the leapfrog scheme.

> stepVelocity :: ( Source a Double
>                 , Source b Double
>                 , Source c Double
>                 )  =>
>                 Array a DIM2 Double ->
>                 Array b DIM2 Double ->
>                 Array c DIM1 Double ->
>                 Array D DIM2 Double
> stepVelocity vs fs ms = Repa.zipWith (+) vs $
>                         Repa.zipWith (/) ms2 $
>                         Repa.zipWith (*) fs dt2
>   where
>     ms2 :: Array D DIM2 Double
>     ms2 = repDim1To2Inner ms
>     dt2 :: Array D DIM2 Double
>     dt2 = extend (Any :. i :. j) $ fromListUnboxed Z [dt]
>     (Z :. i :. j) = extent fs

> stepPosition :: ( Source a Double
>                 , Source b Double
>                 ) =>
>                 Array a DIM2 Double ->
>                 Array b DIM2 Double ->
>                 Array D DIM2 Double
> stepPosition xs vs = Repa.zipWith (+) xs vs

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
> jupiterEccentrity = 0.048775
>
> jupiterMajRad :: Distance
> jupiterMajRad = (jupiterPerihelion + jupiterAphelion) / 2

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

> n :: Double
> n = sqrt $ gConst * sunMass / jupiterMajRad^3

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

> thetaDotP :: Double -- radians per second
> thetaDotP = n * jupiterMajRad^2 * sqrt (1 - jupiterEccentrity^2) / jupiterPerihelion^2
> deltaThetaP :: Double -- radians
> deltaThetaP = thetaDotP * dt / 2
>
> jupiterVPeri :: Speed
> jupiterVPeri = thetaDotP * jupiterPerihelion
>
> jupiterInitX :: Speed
> jupiterInitX = negate $ jupiterVPeri * deltaThetaP
>
> jupiterInitY :: Speed
> jupiterInitY = negate $ jupiterVPeri * (1 - deltaThetaP^2 / 2)
>
> jupiterV :: (Speed, Speed, Speed)
> jupiterV = (jupiterInitX, jupiterInitY, 0.0)

> testParticles :: Array U DIM2 Double
> testParticles = fromListUnboxed (Z :. (4 ::Int) :. spaceDim) [1..12]

> testParticles2 :: Array U DIM2 Double
> testParticles2 = fromListUnboxed (Z :. (2 ::Int) :. spaceDim) [1,1,1,2,2,2]

> testParticles3 :: Array U DIM2 Double
> testParticles3 = fromListUnboxed (Z :. (3 ::Int) :. spaceDim) [1,2,3,5,7,11,13,17,19]

> type Result = IO (Array U DIM3 Double)

> main :: IO ()
> main = do mn <- computeP m :: IO (Array U DIM1 Double)
>           putStrLn $ "Base data..."
>           putStrLn $ show mn
>           putStrLn $ show testParticles3
>           putStrLn $ "Forces..."
>           fsm <- computeP $ forces testParticles3 m :: Result
>           putStrLn $ show fsm
