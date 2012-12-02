{-# OPTIONS_GHC -Wall -fno-warn-name-shadowing -fno-warn-type-defaults #-}

import qualified Data.Vector as V

import Text.PrettyPrint

type Distance = Double
type Mass     = Double
type Force    = Double
type Velocity = Double
type Time     = Double

gConst :: Double
gConst   = 6.674e-11

mSun, mJupiter, mEarth :: Mass
mSun     = 1.9889e30
mJupiter = 1.8986e27
mEarth   = 5.9722e24

dEarth, dSun, dJupiter :: Distance
dEarth   = 1.4960e11
dSun     = 0.0000e11
dJupiter = 7.7855e11

secondsPerDay :: Double
secondsPerDay = 24*60*60

dt :: Time
dt = 10 * secondsPerDay

type PositionV = V.Vector Distance
type VelocityV = V.Vector Velocity

zeroV :: V.Vector Double
zeroV = V.fromList [0.0, 0.0, 0.0]

forcesV :: Int -> V.Vector (V.Vector Distance) -> V.Vector Mass -> V.Vector (V.Vector Force)
forcesV ix positions masses = fs
  where
    v = positions V.! ix
    m = masses    V.! ix

    -- Calculate the pointwise distances between the i-th particle
    -- and the rest of the particles
    ds = V.zipWith (V.zipWith (-)) positions (V.replicate (V.length positions) v)

    -- Calculate the sum of the squares of the distances between the
    -- i-th particle and the rest of the particles and "normalize"
    dsqs = V.map f ds
           where
             f x = (**(3/2)) $ V.sum $ V.map (^2) x

    -- Calculate the forces on the i-th particle
    fs = V.generate (V.length ds) f
         where
           f :: Int -> V.Vector Double
           f i | i == ix   = V.fromList $ take 3 $ repeat 0.0
               | otherwise = V.map g (ds V.! i)
               where
                 g x = gConst * x * (m * (masses V.! i) / (dsqs V.! i))

data Particle = Particle { position :: (Double, Double, Double), mass :: Double }
                  deriving Show

particles :: [Particle]
particles = [ Particle { position = ( 0.0, 0.0, 0.0) , mass = mEarth   }
            , Particle { position = (10.0, 0.0, 0.0) , mass = mSun     }
            , Particle { position = (-10.0, 0.0, 0.0), mass = mJupiter }
            ]

positionVs :: V.Vector (V.Vector Double)
positionVs = V.fromList $ map V.fromList $ map ((\(x, y, z) -> [x, y, z]) . position) particles

massVs :: V.Vector Double
massVs = V.fromList $ map mass particles

data PointedArrayV a = PointedArrayV Int (V.Vector a)
                       deriving Show

instance Functor PointedArrayV where
  fmap f (PointedArrayV i v) = PointedArrayV i (fmap f v)

instance Comonad PointedArrayV where
  coreturn (PointedArrayV i v) = v V.! i
  (PointedArrayV i v) =>> f =
    PointedArrayV i (V.map (f . flip PointedArrayV v) (V.generate (V.length v) id))

f :: PointedArrayV (PositionV, VelocityV) -> (PositionV, VelocityV)
f (PointedArrayV i z) = (newPos, newVel)
                        where
                          x = V.map fst z
                          v = V.map snd z
                          fs = forcesV i x massVs
                          totalForce = V.foldr (V.zipWith (+)) zeroV fs
                          newVel = V.zipWith g totalForce (v V.! i)
                          newPos = V.zipWith h (x V.! i) newVel
                          g force vel  = vel + force * dt / (massVs V.! i)
                          h pos vel    = pos + vel * dt

initPos, initVel :: [V.Vector Double]
initPos =  map V.fromList [[1.496e11, 0, 0], [0, 0, 0]]
initVel =  map V.fromList [[2557.5, 29668.52, 0], [0, 0, 0]]
initPosAndVel :: V.Vector (PositionV, VelocityV)
initPosAndVel =  V.fromList $ zip initPos initVel

prettyPointedArrayV :: PointedArrayV (PositionV, VelocityV) -> String
prettyPointedArrayV (PointedArrayV i x) = render d
    where
      (u, v) = x V.! i
      d = text (show $ V.toList u) <+> text (show $ V.toList v)

testV0, testV1, testV2 :: V.Vector (V.Vector Double)
testV0 = forcesV 0 positionVs massVs
testV1 = forcesV 1 positionVs massVs
testV2 = forcesV 2 positionVs massVs

class Comonad c where
  coreturn :: c a -> a
  (=>>) :: c a -> (c a -> b) -> c b
