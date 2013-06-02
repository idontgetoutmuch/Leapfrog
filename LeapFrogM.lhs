We wish to model planetary motion using the
[leapfrog][leapfrog] method. This is preferred to other numerical methods as it maintains the total energy of the system.

  [leapfrog]: http://en.wikipedia.org/wiki/Leapfrog_integration

In essence, we have to update the position and velocity of each planet half a time step out of phase (hence the name leapfrog) as shown below. 

\begin{align*}
x_i &= x_{i-1} + v_{i - \frac{1/2}}\Delta t \\
a_i &= F(x_i) \\
v_{i + \frac{1}{2]} = v_{i - \frac{1}{2}) + a_i\Delta t
\end{align*}

First some necessary preamble and imports.

> {-# LANGUAGE NoMonomorphismRestriction #-}
> {-# OPTIONS_GHC -Wall -fno-warn-name-shadowing -fno-warn-type-defaults #-}

> import Text.Printf
> import qualified Data.Vector as V

It's easy to confuse units so we define some type synomyms to help. Of
course having these checked by machine would be better.

> type Distance = Double
> type Mass     = Double
> type Force    = Double
> type Speed    = Double
> type Time     = Double

We model the positions and velocities of the planets as vectors of
precisley 3 elements each. Again we would like this invariant to be
checked by machine.

> type PositionV = V.Vector Distance
> type VelocityV = V.Vector Speed

It's not worth giving a type synonym for the gravitational constant
type.

> gConst :: Double -- N (m / kg)^2
> gConst   = 6.674e-11

> secondsPerDay :: Double
> secondsPerDay = 24*60*60

> mSun, mJupiter, mEarth :: Mass
> mSun     = 1.9889e30
> mJupiter = 1.8986e27
> mEarth   = 5.9722e24

The perihelia of the planets in which we are interested.

> dEarthPeri, dSun, dJupiterPeri :: Distance
> dEarthPeri   = 1.496000e11
> dSun         = 0.000000e11
> dJupiterPeri = 7.405736e11

Various other planetary observations.

> dJupiterAp :: Distance
> dJupiterAp        = 8.165208e11
> jupiterEccentrity :: Double
> jupiterEccentrity = 0.048775
> jupiterA :: Distance
> jupiterA          = (dJupiterAp + dJupiterPeri) / 2

FIXME: Explain n!

> n :: Double
> n = sqrt $ gConst * mSun / jupiterA^3

Finally we can calculate Jupiter's velocity by assuming that its
perihelion is on the $x$-axis and that its velocity in the $x$
direction must be $0$.

FIXME: Actually we need the velocity to be one time step before the
perihelion.

> thetaDotP :: Double
> thetaDotP = n * jupiterA^2 * sqrt (1 - jupiterEccentrity^2) / dJupiterPeri^2 -- radians per second
> jupiterV :: (Speed, Speed, Speed)
> jupiterV = (0.0, thetaDotP * dJupiterPeri, 0.0)

The zero vector.

> zeroV :: V.Vector Double
> zeroV = V.fromList [0.0, 0.0, 0.0]

For convenience we collect all the relevant data about a planet into a
record.

> data Planet = Planet { position :: (Distance, Distance, Distance)
>                      , velocity :: (Speed, Speed, Speed)
>                      , mass :: Double
>                      }
>                   deriving Show
>
> planets :: [Planet]
> planets = [ Planet { position = (dEarthPeri, 0.0, 0.0)
>                    , velocity = (2557.5, 29668.52, 0)
>                    , mass = mEarth
>                    }
>           , Planet { position = (dSun,      0.0, 0.0)
>                    , velocity = (0.0, 0.0, 0.0)
>                    , mass = mSun
>                    }
>           , Planet { position = (-dJupiterPeri, 0.0, 0.0)
>                    , velocity = jupiterV
>                    , mass = mJupiter
>                    }
>           ]


We want to calculate with vectors not lists.

> tripleToVec :: (Double, Double, Double) -> V.Vector Double
> tripleToVec (x, y, z) = V.fromList [x, y, z]

> initPos :: [PositionV]
> initPos = take 2 $ map (tripleToVec . position) planets
> initVel :: [VelocityV]
> initVel = take 2 $ map (tripleToVec . velocity) planets

> initPosAndVel :: V.Vector (PositionV, VelocityV)
> initPosAndVel =  V.fromList $ zip initPos initVel

> massVs :: V.Vector Mass
> massVs = V.fromList $ map mass planets

Calculate the force on the $i$-th planet. We could improve the efficiency of this by noting that the force on the $i$-th planet by the $j$-th planet is the same as the force of the $j$-th planet on the $i$-th planet.

> 
> forcesV :: Int -> V.Vector (V.Vector Distance) -> V.Vector Mass -> V.Vector (V.Vector Force)
> forcesV ix positions masses = fs
>   where
>     v = positions V.! ix
>     m = masses    V.! ix

Calculate the pointwise distances between the $i$-th particle
and the rest of the planets. It's a bit of a nuisance to have to explicitly lift arithmetic operators explicitly to operate on vectors; we could avoid this by using typeclasses.

>     ds = V.zipWith (V.zipWith (-)) positions (V.replicate (V.length positions) v)

Calculate the sum of the squares of the distances between the
$i$-th particle and the rest of the planets and "normalize".

>     dsqs = V.map f ds
>            where
>              f x =  (V.sum $ V.map (^2) x) ** (3/2)
> 

Calculate the forces on the $i$-th particle.

>     fs = V.generate (V.length ds) f
>          where
>            f :: Int -> V.Vector Double
>            f i | i == ix   = zeroV
>                | otherwise = V.map g (ds V.! i)
>                where
>                  g x = gConst * x * (m * (masses V.! i) / (dsqs V.! i))

We make the a timeslice of our calculation an array of the state (positions and velocities) of the planets with a distinguished planet.

> data PointedArrayV a = PointedArrayV Int (V.Vector a)
>                        deriving Show

And declare this as a functor and a comonad with the usual comonad defintion.

> instance Functor PointedArrayV where
>   fmap f (PointedArrayV i v) = PointedArrayV i (fmap f v)

> class Comonad c where
>   coreturn :: c a -> a
>   (=>>) :: c a -> (c a -> b) -> c b

> instance Comonad PointedArrayV where
>   coreturn (PointedArrayV i v) = v V.! i
>   (PointedArrayV i v) =>> f =
>     PointedArrayV i (V.map (f . flip PointedArrayV v) (V.generate (V.length v) id))

Now we are ready to move the planets forward in the orbits one step at a time.

> timeStepDays :: Time
> timeStepDays = 10.0

> nTimeSteps :: Int
> nTimeSteps = floor $ 365 / timeStepDays

> dt :: Time
> dt = timeStepDays * secondsPerDay

> oneStep :: PointedArrayV (PositionV, VelocityV) -> (PositionV, VelocityV)
> oneStep (PointedArrayV i z) =
>   (newPos, newVel)
>       where
>         oldPos = V.map fst z
>         oldVel = V.map snd z
> 
>         fs         = forcesV i oldPos massVs
>         totalForce = V.foldr (V.zipWith (+)) zeroV fs
> 
>         newVel = V.zipWith updateVel totalForce (oldVel V.! i)
>         newPos = V.zipWith updatePos (oldPos V.! i) newVel
> 
>         updateVel force vel  = vel + force * dt / (massVs V.! i)
>         updatePos pos vel    = pos + vel * dt

Finally we can run our planetary system.

> test :: IO ()
> test = mapM_ putStrLn $
>        map (f . fst . coreturn) $
>        take nTimeSteps $
>        iterate (=>> oneStep) (PointedArrayV 0 initPosAndVel)
>   where
>     f x = printf "%.5e, %.5e, %.5e" (x V.! 0) (x V.! 1) (x V.! 2)

We can check our solution is working as expected by calculating the
total energy and checking that it remains roughly constant.
