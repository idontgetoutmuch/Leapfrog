{-# OPTIONS_GHC -Wall                     #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing  #-}
{-# OPTIONS_GHC -fno-warn-type-defaults   #-}

{-# LANGUAGE NoMonomorphismRestriction    #-}
{-# LANGUAGE FlexibleContexts             #-}
{-# LANGUAGE ScopedTypeVariables          #-}

module Main (main) where

import qualified Initial as I

import Data.Array.Repa hiding ((++), zipWith)
import qualified Data.Array.Repa as Repa
import Data.Array.Repa.Algorithms.Matrix
import Control.Monad
import Control.Monad.Identity
import Text.Printf
import qualified Data.List as L

stepPositionP :: forall a b c m . ( Monad m
                 , Source a Double
                 , Source b Double
                 , Source c Double
                 ) =>
                 Double ->
                 Array a DIM2 Double ->
                 Array b DIM1 Double ->
                 Array c DIM2 Double ->
                 m (Array U DIM2 Double)
stepPositionP h qs ms ps = do
  do newQs <- computeP $ qs +^ (ps *^ h2 /^ ms2)
     return newQs
    where
      (Z :. i :. j) = extent ps

      h2  = extend (Any :. i :. j) $ fromListUnboxed Z [h]
      ms2 = extend (Any :. j) ms

stepMomentumP :: forall a b c m . ( Monad m
                 , Source a Double
                 , Source b Double
                 , Source c Double
                 ) =>
                 Double ->
                 Double ->
                 Array a DIM2 Double ->
                 Array b DIM1 Double ->
                 Array c DIM2 Double ->
                 m (Array U DIM2 Double)
stepMomentumP gConst h qs ms ps =
  do fs <- sumP $ transpose $ zeroDiags' fss
     newPs <- computeP $ ps +^ (fs *^ dt2)
     return newPs
  where
    is = repDim2to3Outer $ prodPairsMasses ms
    qDiffs = pointDiffs qs
    preDs = Repa.map (^3) $
            Repa.map sqrt $
            sumS $
            Repa.map (^2) $
            qDiffs
    ds    = repDim2to3Outer preDs
    preFs = Repa.map (* (negate gConst)) $
            qDiffs /^ ds
    fss = is *^ preFs
    Z :.i :. j :. k = extent fss
    dt2             = extend (Any :. i :. k) $ fromListUnboxed Z [h]

stepOnceP :: ( Monad m
             , Source a Double
             , Source b Double
             , Source c Double
             ) =>
             Double ->
             Double ->
             Array b DIM1 Double ->
             Array a DIM2 Double ->
             Array c DIM2 Double ->
             m (Array U DIM2 Double, Array U DIM2 Double)
stepOnceP gConst h ms qs ps = do
  newPs <- stepMomentumP gConst h qs ms ps
  newQs <- stepPositionP h qs ms newPs
  return (newQs, newPs)

-- FIXME
repDim2to3Outer a = extend (Any :. I.spaceDim) a


type Momenta   = Array U DIM2 Double
type Positions = Array U DIM2 Double
type Masses    = Array U DIM1 Double

zeroDiags x = traverse x id f
  where
    f _ (Z :. i :. j) | i == j    = 0.0
                      | otherwise = x!(Z :. i :. j)
                                    
zeroDiags' x = traverse x id f
  where
    f _ (Z :. i :. j :. k) | i == j    = 0.0
                           | otherwise = x!(Z :. i :. j :. k)

hamiltonianP :: Double -> Masses -> Momenta -> Positions -> IO Double
hamiltonianP gConst ms qs ps = do
  preKes <- sumP $ ps *^ ps
  ke     <- sumP $ preKes /^ ms

  ds2 <- sumP $ Repa.map (^2) $ pointDiffs qs
  let ds   = Repa.map sqrt ds2
      is   = prodPairsMasses ms
      pess = zeroDiags $ Repa.map (* (negate gConst)) $ is /^ ds
  pes <- sumP pess
  pe  <- sumP pes
  te :: Array U DIM0 Double <- computeP $ ke +^ pe
  return $ head $ toList $ Repa.map (* 0.5) te

repDim1To2Outer :: Source a Double =>
                   Array a DIM1 Double ->
                   Array D DIM2 Double
repDim1To2Outer a = extend (Any :. i :. All) a
  where (Z :. i) = extent a

prodPairsMasses :: Source a Double =>
                   Array a DIM1 Double ->
                   Array D DIM2 Double
prodPairsMasses ms = ns *^ (transpose ns)
  where
    ns = repDim1To2Outer ms

transposeOuter :: Source a Double =>
              Array a DIM3 Double ->
              Array D DIM3 Double
transposeOuter qs = backpermute (f e) f qs
  where
    e = extent qs
    f (Z :. i :. i' :. j) = Z :. i' :. i :. j

pointDiffs :: Source a Double =>
              Array a DIM2 Double ->
              Array D DIM3 Double
pointDiffs qs = qss -^ (transposeOuter qss)
  where qss = replicateRows qs

replicateRows :: Source a Double =>
                 Array a DIM2 Double ->
                 Array D DIM3 Double
replicateRows a = extend (Any :. i :. All) a
  where (Z :. i :. _j) = extent a

type Distance = Double
type Mass     = Double
type Force    = Double
type Speed    = Double
type Energy   = Double

initPs :: Array U DIM2 Double
initPs = runIdentity $ computeP $ ms2 *^ initVs
  where
    (Z :. i :. j) = extent initVs
    ms2 = extend (Any :. j) masses

initRs :: Array U DIM2 Distance
initRs = fromListUnboxed (Z :. nBodies :. spaceDim) $ concat xs
  where
    nBodies = length xs
    xs = [ [earthX,   earthY,   earthZ]
         , [jupiterX, jupiterY, jupiterZ]
         , [sunX,     sunY,     sunZ]
         ]
    (earthX,   earthY,   earthZ)   = earthR
    (jupiterX, jupiterY, jupiterZ) = jupiterR
    (sunX,     sunY,     sunZ)     = sunR

masses :: Array U DIM1 Mass
masses = fromListUnboxed (Z :. nBodies) xs
  where
    nBodies = length xs
    xs = [ earthMass
         , jupiterMass
         , sunMass
         ]

stepN :: forall m . Monad m =>
         Int -> Double -> Double -> Masses -> Positions -> Momenta ->
         m (Positions, Momenta)
stepN n gConst dt masses = curry updaterMulti
  where
    updaterMulti = foldr (>=>) return updaters
    updaters :: [(Positions, Momenta) -> m (Positions, Momenta)]
    updaters = replicate n (uncurry (stepOnceP gConst dt masses))



