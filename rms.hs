{-# LANGUAGE NoMonomorphismRestriction #-}

import Data.Array.Repa as Repa

import Data.Random ()
import Data.Random.Distribution.Uniform
import Data.RVar

import System.Random
import Control.Monad
import Control.Monad.State

space :: Int
space = 3

g = 6.67384e-11             -- gravitational constant
max_tolerance = 0.1         -- allow 10% change in energy
k = 24*60*60                -- seconds in a day
   
-- Initial values
r = 1e1                    -- radius of sphere particles contained in
eps = 0.1*r                 -- softening constant
days = 36500*100            -- total time in days
t = days*k                  -- total time 
timestep_days = 10          -- timestep in days
dt = timestep_days*k        -- timestep 
nIters = floor (t/dt)       -- number of iterations

nBodies :: Int
nBodies = 4                -- number of bodies 
mass = 1e24                 -- mass

-- array of equal masses
m = Repa.map (mass*) $ 
    fromListUnboxed (Z :. nBodies) $
    take nBodies $ repeat 1.0

rands :: Int -> Double -> Double -> [Double]
rands n a b =
  fst $ runState (replicateM n (sampleRVar (uniform a b))) (mkStdGen seed)
    where
      seed = 0

-- The random position is in a box r^3
randomPoss :: Int -> Double -> Array U (Z :. Int :. Int) Double
randomPoss n r = fromListUnboxed (Z :. n :. space) $
                 Prelude.map (fromIntegral . round) $
                 rands (3*n) 0 r 

startPoss = randomPoss nBodies r

distances ps = computeP $ Repa.map sqrt $ sumS $ Repa.map (^2) ps

forces ps ms = undefined
  where
    ds = Repa.zipWith (-) ps (twizzle ps)

twizzle ps = traverse ps id (\f (Z :. i :. j) -> f (Z :. (i + 1) `mod` nBodies :. j))

positions = Repa.zipWith (-) startPoss startPoss

cube :: Array D (Z :.Int :. Int :. Int) Double
cube = extend (Any :. i :. All) startPoss
         where (Z :. i :. j) = extent startPoss

cube' :: Array U (Z :. Int :. Int) Double -> Array D (Z :.Int :. Int :. Int) Double
cube' a = extend (Any :. i :. All) a
            where (Z :. i :. j) = extent a

twizzleds ps = backpermute e f ps
                 where
                   e = extent ps
                   f (Z :. i :. i' :. j) = Z :. (i + 1) `mod` nBodies :. i' :. j

main = do mn <- computeP m :: IO (Array U (Z :. Int) Double)
          putStrLn $ show mn
          dm <- distances startPoss :: IO (Array U (Z :. Int) Double)
          putStrLn $ show dm
          pm <- computeP positions :: IO (Array U (Z :. Int :. Int) Double)
          putStrLn $ show pm
          tm <- computeP $ twizzle startPoss :: IO (Array U (Z :. Int :. Int) Double)
          putStrLn $ show tm
          putStrLn $ show startPoss
          em <- computeP $ extend (Any :. (5 :: Int) :. All) (randomPoss 3 10) :: IO (Array U (Z :. Int :. Int :. Int) Double)
          putStrLn $ show em
          cm <- computeP cube :: IO (Array U (Z :. Int :. Int :. Int) Double)
          putStrLn $ show cm

foo :: IO (Array U (Z :. Int) Double)
foo = computeP $ slice (randomPoss nBodies r) (Z :. (0 :: Int) :. All)

bar ps = traverse ps (\x@(Z :. i :. j) -> (Z :. i))
                     (\f (Z :. i) -> f (Z :.i :. (0 :: Int)))

baz :: IO (Array U (Z :. Int) Double)
baz = computeP $ bar startPoss