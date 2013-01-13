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

extraDim :: Array U (Z :. Int :. Int) Double -> Array D (Z :.Int :. Int :. Int) Double
extraDim a = extend (Any :. i :. All) a
            where (Z :. i :. j) = extent a

extraDim' :: Array D (Z :. Int :. Int) Double -> Array D (Z :.Int :. Int :. Int) Double
extraDim' a = extend (Any :. (3 :: Int)) a
            where (Z :. i :. j) = extent a

transpose3 ps = backpermute (f e) f ps
                   where
                     e = extent ps
                     f (Z :. i :. i' :. j) = Z :. i' :. i :. j

pointDiffs ps = Repa.zipWith (-) qs (transpose3 qs)
                  where qs = extraDim ps

distances = f
  where f = Repa.map sqrt . foldS (+) 0 . Repa.map (^2) . pointDiffs

forces qs = Repa.map (* (negate g)) $ Repa.zipWith (/) (pointDiffs qs) ds
  where ds = extraDim' $
             Repa.map sqrt $
             Repa.map (+ (eps^2)) $
             foldS (+) 0 $
             Repa.map (^2) $
             pointDiffs qs

testParticles :: Array U (Z :. Int :. Int) Double
testParticles = fromListUnboxed (Z :. (4 ::Int) :. (3 :: Int)) [1..12]

testParticles2 :: Array U (Z :. Int :. Int) Double
testParticles2 = fromListUnboxed (Z :. (2 ::Int) :. (3 :: Int)) [1,2,3,5,7,11]

testParticles3 :: Array U (Z :. Int :. Int) Double
testParticles3 = fromListUnboxed (Z :. (3 ::Int) :. (3 :: Int)) [1,2,3,5,7,11,13,17,19]

type Result = IO (Array U (Z :. Int :. Int :. Int) Double)

main = do mn <- computeP m :: IO (Array U (Z :. Int) Double)
          putStrLn $ show mn
          -- dm <- distances startPoss :: IO (Array U (Z :. Int) Double)
          -- putStrLn $ show dm
          -- pm <- computeP positions :: IO (Array U (Z :. Int :. Int) Double)
          -- putStrLn $ show pm

          putStrLn $ show testParticles2
          tpcm <- computeP $ extraDim testParticles2 :: Result
          putStrLn $ show tpcm
          rcm <- computeP $ transpose3 $ extraDim $ testParticles2 :: Result
          putStrLn $ show rcm
          diffs <- computeP $ pointDiffs testParticles2 :: Result
          putStrLn $ show diffs
          dsm <- computeP $ distances testParticles2 :: IO (Array U (Z :. Int :. Int) Double)
          putStrLn $ show dsm
          fsm <- computeP $ forces testParticles2 :: Result
          putStrLn $ show fsm
