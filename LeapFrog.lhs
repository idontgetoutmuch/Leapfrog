First some necessary preamble and imports

> {-# LANGUAGE NoMonomorphismRestriction #-}
> {-# OPTIONS_GHC -Wall -fno-warn-name-shadowing -fno-warn-type-defaults #-}

> import Text.Printf
> import qualified Data.Vector as V

It's easy to confuse units so we define some type synomyms to help (of
course having these checked by machine would be better).

> type Distance = Double
> type Mass     = Double
> type Force    = Double
> type Velocity = Double
> type Time     = Double

We model the positions and velocities of the planets as vectors of
precisley 3 elements each. Again we would like this invariant to be
checked by machine.

> type PositionV = V.Vector Distance
> type VelocityV = V.Vector Velocity

No type synonym for the gravitational constant type.

> gConst :: Double -- N (m / kg)^2
> gConst   = 6.674e-11

> secondsPerDay :: Double
> secondsPerDay = 24*60*60

> mSun, mJupiter, mEarth :: Mass
> mSun     = 1.9889e30
> mJupiter = 1.8986e27
> mEarth   = 5.9722e24

The perihelia of the planets in which we are interested.

> dEarth, dSun, dJupiter :: Distance
> dEarth   = 1.4960e11
> dSun     = 0.0000e11
> dJupiter = 7.7855e11

> zeroV :: V.Vector Double
> zeroV = V.fromList [0.0, 0.0, 0.0]

> initPos :: [PositionV]
> initPos = take 2 $ positionVs
> initVel :: [VelocityV]
> initVel =  map V.fromList [[2557.5, 29668.52, 0], [0, 0, 0]]
> initPosAndVel :: V.Vector (PositionV, VelocityV)
> initPosAndVel =  V.fromList $ zip initPos initVel

> data Particle = Particle { position :: (Double, Double, Double)
>                          , mass :: Double
>                          }
>                   deriving Show

> particles :: [Particle]
> particles = [ Particle { position = (dEarth,    0.0, 0.0), mass = mEarth   }
>             , Particle { position = (dSun,      0.0, 0.0), mass = mSun     }
>             , Particle { position = (-dJupiter, 0.0, 0.0), mass = mJupiter }
>             ]

> positionVs :: [V.Vector Double]
> positionVs = map V.fromList $ map ((\(x, y, z) -> [x, y, z]) . position) particles
> 
> massVs :: V.Vector Double
> massVs = V.fromList $ map mass particles
> 
> forcesV :: Int -> V.Vector (V.Vector Distance) -> V.Vector Mass -> V.Vector (V.Vector Force)
> forcesV ix positions masses = fs
>   where
>     v = positions V.! ix
>     m = masses    V.! ix
> 
>     -- Calculate the pointwise distances between the i-th particle
>     -- and the rest of the particles
>     ds = V.zipWith (V.zipWith (-)) positions (V.replicate (V.length positions) v)
> 
>     -- Calculate the sum of the squares of the distances between the
>     -- i-th particle and the rest of the particles and "normalize"
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

test :: IO ()
test = mapM_ putStrLn $
       map (f . fst . coreturn) $
       take nTimeSteps $
       iterate (=>> oneStep) (PointedArrayV 0 initPosAndVel)
  where
    f x = printf "%.5e, %.5e, %.5e" (x V.! 0) (x V.! 1) (x V.! 2)
\end{code}