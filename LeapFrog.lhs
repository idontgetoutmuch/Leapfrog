First some necessary preamble and imports

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

No type synonym for the gravitational constant type

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

> dJupiterAp        = 8.165208e11
> jupiterEccentrity = 0.048775
> jupiterA          = (dJupiterAp + dJupiterPeri) / 2

FIXME: Explain n!

> n :: Double
> n = sqrt $ gConst * mSun / jupiterA^3

Finally we can calculate Jupiter's velocity by assuming that its perihelion on the $x$-axis and then its velocity in the $x$ direction must be $0$.

> thetaDotP = n * jupiterA^2 * sqrt (1 - jupiterEccentrity^2) / dJupiterPeri^2 -- radians per second
> jupiterV = (0.0, thetaDotP * dJupiterPeri, 0.0)

The zero vector.

> zeroV :: V.Vector Double
> zeroV = V.fromList [0.0, 0.0, 0.0]

> data Planet = Planet { position :: (Double, Double, Double)
>                      , velocity :: (Speed, Speed, Speed)
>                      , mass :: Double
>                      }
>                   deriving Show

> planets :: [Planet]
> planets = [ Planet { position = (dEarthPeri,    0.0, 0.0)
>                    , velocity = (2557.5, 29668.52, 0)
>                    , mass = mEarth   }
>           , Planet { position = (dSun,      0.0, 0.0)
>                    , velocity = (0.0, 0.0, 0.0)
>                    , mass = mSun     }
>           , Planet { position = (-dJupiterPeri, 0.0, 0.0), mass = mJupiter }
>           ]

> initPos :: [PositionV]
> initPos = take 2 $ positionVs
> initVel :: [VelocityV]
> initVel =  map V.fromList [[2557.5, 29668.52, 0], [0, 0, 0]]
> initPosAndVel :: V.Vector (PositionV, VelocityV)
> initPosAndVel =  V.fromList $ zip initPos initVel

> positionVs :: [V.Vector Double]
> positionVs = map V.fromList $ map ((\(x, y, z) -> [x, y, z]) . position) planets
> 
> massVs :: V.Vector Double
> massVs = V.fromList $ map mass planets
> 
> forcesV :: Int -> V.Vector (V.Vector Distance) -> V.Vector Mass -> V.Vector (V.Vector Force)
> forcesV ix positions masses = fs
>   where
>     v = positions V.! ix
>     m = masses    V.! ix
> 
>     -- Calculate the pointwise distances between the i-th particle
>     -- and the rest of the planets
>     ds = V.zipWith (V.zipWith (-)) positions (V.replicate (V.length positions) v)
> 
>     -- Calculate the sum of the squares of the distances between the
>     -- i-th particle and the rest of the planets and "normalize"
>     dsqs = V.map f ds
>            where
>              f x = (**(3/2)) $ V.sum $ V.map (^2) x
> 
>     -- Calculate the forces on the i-th particle
>     fs = V.generate (V.length ds) f
>          where
>            f :: Int -> V.Vector Double
>            f i | i == ix   = zeroV
>                | otherwise = V.map g (ds V.! i)
>                where
>                  g x = gConst * x * (m * (masses V.! i) / (dsqs V.! i))

> data PointedArrayV a = PointedArrayV Int (V.Vector a)
>                        deriving Show

> instance Functor PointedArrayV where
>   fmap f (PointedArrayV i v) = PointedArrayV i (fmap f v)

> class Comonad c where
>   coreturn :: c a -> a
>   (=>>) :: c a -> (c a -> b) -> c b

> instance Comonad PointedArrayV where
>   coreturn (PointedArrayV i v) = v V.! i
>   (PointedArrayV i v) =>> f =
>     PointedArrayV i (V.map (f . flip PointedArrayV v) (V.generate (V.length v) id))

> timeStepDays :: Time
> timeStepDays = 0.01

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

> test :: IO ()
> test = mapM_ putStrLn $
>        map (f . fst . coreturn) $
>        take nTimeSteps $
>        iterate (=>> oneStep) (PointedArrayV 0 initPosAndVel)
>   where
>     f x = printf "%.5e, %.5e, %.5e" (x V.! 0) (x V.! 1) (x V.! 2)

\end{code}