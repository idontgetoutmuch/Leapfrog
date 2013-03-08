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

> import Debug.Trace

We work in a 3 dimensional Euclidean space.

> spaceDim :: Int
> spaceDim = 3

Some physical constants for our system.

> g, k :: Double
> g = 6.67384e-11             -- gravitational constant
> k = 24*60*60                -- seconds in a day
> r = 1e1                    -- radius of sphere particles contained in

> nBodies :: Int
> nBodies = 3                -- number of bodies
> mass = 1e24                 -- mass

And some constants for our finite difference method.

> eps = 0.1*r                 -- softening constant
> days = 36500*100            -- total time in days
> t = days*k                  -- total time
> timestep_days = 10          -- timestep in days
> dt = timestep_days*k        -- timestep
> nIters = floor (t/dt)       -- number of iterations

We represent our system of bodies of identical mass as a
1-dimensional (unboxed) array.

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
> startPoss = randomPoss nBodies r

We wish to calculate a two dimensional array of forces where each
force is itself a one dimensional vector of three elements (in the
case of the three dimensional Euclidean space in which we are
operating).

$$
F_i = G\sum_{j\neq i} -m_i m_j\frac{\vec{r}_i - \vec{r}_j}{\abs{vec{r}_i - vec{r}_j}^3}
$$

> replicateRows :: Source a Double =>
>             Array a DIM2 Double -> Array D DIM3 Double
> replicateRows a = extend (Any :. i :. All) a
>             where (Z :. i :. j) = extent a

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



> extraDim' :: Source a Double =>
>              Array a DIM2 Double -> Array D DIM3 Double
> extraDim' a = extend (Any :. spaceDim) a
>
> transposeInner :: Source a Double =>
>               Array a DIM3 Double -> Array D DIM3 Double
> transposeInner ps = backpermute (f e) f ps
>                       where
>                         e = extent ps
>                         f (Z :. i :. i' :. j) = Z :. i' :. i :. j

> pointDiffs :: Source a Double =>
>               Array a DIM2 Double -> Array D DIM3 Double
> pointDiffs ps = Repa.zipWith (-) qs (transposeInner qs)
>                   where qs = replicateRows ps
>
> forces :: (Source a Double, Source b Double) =>
>            Array a DIM2 Double -> Array b DIM1 Double ->
>            Array D DIM3 Double
> forces qs ms = Repa.zipWith (*) fs is
>   where ds = extraDim' $
>              Repa.map sqrt $
>              Repa.map (+ (eps^2)) $
>              foldS (+) 0 $
>              Repa.map (^2) $
>              pointDiffs qs
>         is = extraDim' $ prodPairsMasses ms
>         fs = Repa.map (* (negate g)) $ Repa.zipWith (/) (pointDiffs qs) ds

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
