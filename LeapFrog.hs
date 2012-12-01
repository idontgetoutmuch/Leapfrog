{-# LANGUAGE NoMonomorphismRestriction, TypeFamilies #-}

import qualified Data.Array.Repa as Repa
import Data.Array.Repa (U(..), Z(..), (:.)(..), (!), All(..), DIM1, DIM2, DIM3, Array)

import qualified Data.Vector as V

import Text.PrettyPrint
import Debug.Trace

type Distance = Double
type Mass     = Double
type Force    = Double
type Velocity = Double

-- Gravitational constant
-- gConst :: Double
gConst   = 6.674e-11
mSun     = 1.9889e30
mJupiter = 1.8986e27
mEarth   = 5.9722e24

dEarth   = 1.4960e11
dSun     = 0.0000e11
dJupiter = 7.7855e11

-- Time step
dt = 0.1

-- forces :: Int -> Array U DIM2 Distance ->
--              Array U DIM1 Mass ->
--              Array Repa.D DIM2 Double
forces ix positions masses = fs
  where
    -- Get the dimensions of the array
    -- dx, dy :: Int
    Z :. dx :. dy = Repa.extent positions
    
    -- Get the i-th row
    -- v :: Array Repa.D DIM1 Double
    v  = Repa.slice positions (Z :. ix :. All)

    -- m :: Double
    m = masses!(Z :. ix)

    -- Remove the i-th row
    -- vs :: Array Repa.D DIM2 Double
    vs = Repa.fromFunction (Z :. dx - 1 :. dy)
                           (\(Z :. jx :. kx) -> positions ! (Z :. f jx :. kx))
           where
             f jx | jx < ix   = jx
                  | otherwise = jx + 1

    -- ms :: Array Repa.D DIM1 Double
    ms = Repa.fromFunction (Z :. dx - 1)
                               (\(Z :. jx) -> masses ! (Z :. f jx))
               where
                 f jx | jx < ix   = jx
                      | otherwise = jx + 1

    -- Calculate the pointwise distances between the i-th particle
    -- and the rest of the particles
    -- ds :: Array Repa.D DIM2 Double
    ds = Repa.traverse vs id (\f (Z :. ix :. jx) -> g jx $ f (Z :. ix :. jx))
           where
             g jx x = x - v!(Z :. jx)

    -- Calculate the sum of the squares of the distances between the
    -- i-th particle and the rest of the particles and "normalize"
    -- dsqs :: Array Repa.D DIM1 Double
    dsqs = Repa.map (**(3/2)) $ Repa.foldS (+) 0 $ Repa.map (^2) ds

    -- Calculate the forces on the i-th particle
    -- fs :: Array Repa.D DIM2 Double
    fs = Repa.traverse ds id (\f (Z :. iy :. jy) -> g iy $ f (Z :. iy :. jy))
           where
             g iy x = gConst * x * m * ms!(Z :. iy) / dsqs!(Z :. iy)

forcesV :: Int -> V.Vector (V.Vector Double) -> V.Vector Double -> V.Vector (V.Vector Double)
forcesV ix positions masses = fs
  where
    -- Get the i-th row
    v :: V.Vector Double
    v = positions V.! ix

    m :: Double
    m = masses V.! ix

    -- Remove the i-th row
    -- vs :: V.Vector (V.Vector Double)
    -- vs = sliceDrop ix positions
    vs :: V.Vector (V.Vector Double)
    vs = positions

    sliceDrop ix ws =
      if V.null ys
        then xs
        else xs V.++ (V.tail ys)
      where
        (xs, ys) = V.splitAt ix ws
        
    ms :: V.Vector Double
    -- ms = sliceDrop ix masses
    ms = masses

    -- Calculate the pointwise distances between the i-th particle
    -- and the rest of the particles
    ds :: V.Vector (V.Vector Double)
    ds = V.zipWith (V.zipWith (-)) vs (V.replicate (V.length vs) v)

    -- Calculate the sum of the squares of the distances between the
    -- i-th particle and the rest of the particles and "normalize"
    dsqs :: V.Vector Double
    dsqs = V.map f ds
           where
             f :: V.Vector Double -> Double
             f x = (**(3/2)) $ V.sum $ V.map (^2) x

    -- Calculate the forces on the i-th particle
    fs :: V.Vector (V.Vector Double)
    fs = V.generate (V.length ds) f
         where
           f :: Int -> V.Vector Double
           f i | i == ix   = V.fromList $ take 3 $ repeat 0.0
               | otherwise = V.map g (ds V.! i)
               where
                 g x = gConst * x * (m * (ms V.! i) / (dsqs V.! i))

data Particle = Particle { position :: (Double, Double, Double), mass :: Double }
                  deriving Show

particles = [ Particle { position = ( 0.0, 0.0, 0.0), mass = 2.0 }
            , Particle { position = (10.0, 0.0, 0.0), mass = 3.0 }
            , Particle { position = (-10.0, 0.0, 0.0), mass = 3.0 }
            ]

positions :: Array U DIM2 Double
positions = Repa.fromListUnboxed (Z :. (3 :: Int) :. (3 :: Int)) $
             concatMap ((\(x, y, z) -> [x, y, z]) . position) particles

positionVs :: V.Vector (V.Vector Double)
positionVs = V.fromList $ map V.fromList $ map ((\(x, y, z) -> [x, y, z]) . position) particles

masses :: Array U DIM1 Double
masses = Repa.fromListUnboxed (Z :. 3) $
         map mass particles

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

safeIndex :: Show a => Int -> V.Vector a -> String -> a
safeIndex i v label = if (i < 0) || (i > V.length v - 1)
                      then
                          error $ label ++ " " ++ (show v) ++ " " ++ show i
                      else
                          v V.! i

f :: PointedArrayV (PositionV, VelocityV) -> (PositionV, VelocityV)
f (PointedArrayV i z) = trace ("Index:\n" ++ show i ++ "\n" ++
                               "New Vel:\n" ++ show w ++ "\n" ++
                               "Forces:\n" ++ show fs) $ (u, w)
                        where
                          x = V.map fst z
                          v = V.map snd z
                          fs :: V.Vector (V.Vector Double)
                          fs = forcesV i x massVs
                          f = V.foldr (V.zipWith (+)) (V.fromList [0.0, 0.0, 0.0]) fs
                          w = V.zipWith g f {- (fs V.! i) -} (v V.! i)
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

test0 = Repa.computeP $ forces 0 positions masses :: IO (Array U DIM2 Double)
testV0 = forcesV 0 positionVs massVs
test1 = Repa.computeP $ forces 1 positions masses :: IO (Array U DIM2 Double)
testV1 = forcesV 1 positionVs massVs
test2 = Repa.computeP $ forces 2 positions masses :: IO (Array U DIM2 Double)
testV2 = forcesV 2 positionVs massVs

-- For each time step we have a 3 by n array. Thus if we have m time
-- steps we have a m by 3 by n array.

vGrid :: Array U DIM1 Mass -> Array U DIM2 Force -> Array U DIM3 Velocity
vGrid m f = undefined

class Foo c where
  type Bar f :: * -> *
  baz :: c a -> Bar c a
  cor :: c a -> (c a -> Bar c b) -> c b

data PointedArray a = PointedArray Int (Array Repa.D DIM2 a)

instance Functor PointedArray where
  fmap f (PointedArray i x) = PointedArray i (Repa.map f x)

instance Foo PointedArray where
  type Bar PointedArray = Array Repa.D DIM1
  baz (PointedArray i x) = Repa.slice x (Z :. i :. All)
  cor (PointedArray i x) f = PointedArray i undefined

x = Repa.fromListUnboxed (Z :. (5 :: Int) :. (3 :: Int)) [1..15]

class Comonad c where
  coreturn :: c a -> a
  (=>>) :: c a -> (c a -> b) -> c b

data PointedArray' a = PointedArray' Int (V.Vector a)
  deriving Show

instance Functor PointedArray' where
  fmap f (PointedArray' i x) = PointedArray' i (fmap f x)

instance Comonad PointedArray' where
  coreturn (PointedArray' i x) = x V.! i
  (PointedArray' i x) =>> f =
    PointedArray' i (V.map (f . flip PointedArray' x) (V.generate (V.length x) id))

-- tc :: PointedArray a -> (PointedArray a -> Array Repa.D DIM1 a) -> PointedArray a
-- tc (PointedArray i x) f =
--   PointedArray i (Repa.map (f . flip PointedArray x) undefined)


xs = Repa.fromListUnboxed (Z :. 3) [1, 2, 3]

removeOne ix xs = Repa.fromFunction (Z :. dx - 1) (\(Z :. jx) -> xs ! (Z :. f jx))
       where
         Z :. dx = Repa.extent xs
         f jx | jx < ix   = jx
              | otherwise = jx + 1

test = Repa.computeP $ removeOne 1 xs :: IO (Array U DIM1 Float)