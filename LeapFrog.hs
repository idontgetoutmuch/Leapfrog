{-# LANGUAGE NoMonomorphismRestriction, TypeFamilies #-}

import qualified Data.Array.Repa as Repa
import Data.Array.Repa (U(..), Z(..), (:.)(..), (!), All(..), DIM1, DIM2, DIM3, Array)

import qualified Data.Vector as V

type Distance = Double
type Mass     = Double
type Force    = Double
type Velocity = Double

-- Gravitational constant
-- gConst :: Double
gConst = {- 1.0 --} 6.674e-11
mSun = 1.9889e30
dEarth = 1.4960e11

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

data Particle = Particle { position :: (Double, Double, Double), mass :: Double }
                  deriving Show

particles = [ Particle { position = ( 0.0, 0.0, 0.0), mass = 2.0 }
            , Particle { position = (10.0, 0.0, 0.0), mass = 3.0 }
            , Particle { position = (-10.0, 0.0, 0.0), mass = 3.0 }
            ]

positions :: Array U DIM2 Double
positions = Repa.fromListUnboxed (Z :. (3 :: Int) :. (3 :: Int)) $
             concatMap ((\(x, y, z) -> [x, y, z]) . position) particles

masses :: Array U DIM1 Double
masses = Repa.fromListUnboxed (Z :. 3) $
         map mass particles

-- test0 = Repa.computeP $ forces 0 positions masses :: IO (Array U DIM2 Double)
-- test1 = Repa.computeP $ forces 1 positions masses :: IO (Array U DIM2 Double)
-- test2 = Repa.computeP $ forces 2 positions masses :: IO (Array U DIM2 Double)

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

tc :: PointedArray a -> (PointedArray a -> Array Repa.D DIM1 a) -> PointedArray a
tc (PointedArray i x) f =
  PointedArray i (Repa.map (f . flip PointedArray x) undefined)