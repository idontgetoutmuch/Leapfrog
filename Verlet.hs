{-# LANGUAGE ScopedTypeVariables, BangPatterns, UnboxedTuples #-}

import qualified Numeric.LinearAlgebra as LA
import qualified Data.List             as L


type Parameters = [Double]
type Density    = Parameters -> Double
type Gradient   = Parameters -> Parameters
type Particle   = (Parameters, Parameters)

leapfrog :: Gradient -> Particle -> Double -> Particle
leapfrog glTarget (t, r) epsilon= (tf, rf)
  where
    rm = adjustMomentum glTarget epsilon t r
    tf = adjustPosition epsilon rm t
    rf = adjustMomentum glTarget epsilon tf rm

adjustMomentum :: Fractional c => (t -> [c]) -> c -> t -> [c] -> [c]
adjustMomentum glTarget e t r = r .+ ((e / 2) .* glTarget t)

adjustPosition :: Num c => c -> [c] -> [c] -> [c]
adjustPosition e r t = t .+ (e .* r)

symplecticEuler :: (LA.Vector Double -> LA.Vector Double) ->
                   LA.Vector Double ->
                   LA.Vector Double ->
                   LA.Matrix Double
symplecticEuler f pms ts = LA.fromRows $
                           map snd $
                           snd selectedTimeSteps
  where
    selectedTimeSteps :: ([(Double, LA.Vector Double)],
                          [(Double, LA.Vector Double)])
    selectedTimeSteps =
      L.mapAccumL (\s t -> let s' = dropWhile (\(x, _) -> x < t) s in (s', head s'))
                  allTimeSteps vs

    allTimeSteps :: [(Double, LA.Vector Double)]
    allTimeSteps = zip us (iterate stepOnce pms)

    us :: [Double]
    us = 0.0 : map (+hi) us
    vs :: [Double]
    vs = LA.toList ts

    stepOnce :: LA.Vector Double ->
                LA.Vector Double
    stepOnce pms = LA.vjoin [newPs, newMs]
      where
        hs = LA.fromList $ replicate (2*n) hi
        deltas = hs * f pms
        newMs = LA.subVector n n pms + LA.subVector n n deltas
        oldPs = LA.subVector 0 n pms
        newDeltas = hs * f (LA.vjoin [oldPs, newMs])
        newPs = LA.subVector 0 n pms + LA.subVector 0 n newDeltas
