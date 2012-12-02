import qualified Data.Vector as V

import Text.PrettyPrint
import Debug.Trace

type Distance = Double
type Mass     = Double
type Force    = Double
type Velocity = Double

gConst   = 6.674e-11
mSun     = 1.9889e30
mJupiter = 1.8986e27
mEarth   = 5.9722e24

dEarth   = 1.4960e11
dSun     = 0.0000e11
dJupiter = 7.7855e11

secondsPerDay = 24*60*60
dt            = 10 * secondsPerDay

forcesV :: Int -> V.Vector (V.Vector Double) -> V.Vector Double -> V.Vector (V.Vector Double)
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

type PositionV = V.Vector Double
type VelocityV = V.Vector Double

f :: PointedArrayV (PositionV, VelocityV) -> (PositionV, VelocityV)
f (PointedArrayV i z) = (u, w)
                        where
                          x = V.map fst z
                          v = V.map snd z
                          fs :: V.Vector (V.Vector Double)
                          fs = forcesV i x massVs
                          f = V.foldr (V.zipWith (+)) (V.fromList [0.0, 0.0, 0.0]) fs
                          w = V.zipWith g f (v V.! i)
                          u = V.zipWith h (x V.! i) w
                          g force vel  = vel + force * dt / (massVs V.! i)
                          h pos vel    = pos + vel * dt

initPos :: [V.Vector Double]
initPos =  map V.fromList [[1.496e11, 0, 0], [0, 0, 0]]

initVel :: [V.Vector Double]
initVel =  map V.fromList [[2557.5, 29668.52, 0], [0, 0, 0]]

initPosAndVel :: V.Vector (PositionV, VelocityV)
initPosAndVel =  V.fromList $ zip initPos initVel

prettyPointedArrayV :: PointedArrayV (PositionV, VelocityV) -> String
prettyPointedArrayV (PointedArrayV i x) = render d
    where
      (u, v) = x V.! i
      d = text (show $ V.toList u) <+> text (show $ V.toList v)

testV0 = forcesV 0 positionVs massVs
testV1 = forcesV 1 positionVs massVs
testV2 = forcesV 2 positionVs massVs

class Comonad c where
  coreturn :: c a -> a
  (=>>) :: c a -> (c a -> b) -> c b
