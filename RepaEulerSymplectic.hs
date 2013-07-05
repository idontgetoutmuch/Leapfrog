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
import Control.Monad

stepPosition :: forall a b c m . ( Monad m
                 , Source a Double
                 , Source b Double
                 , Source c Double
                 ) =>
                 Double ->
                 Array a DIM2 Double ->
                 Array b DIM1 Double ->
                 Array c DIM2 Double ->
                 m (Array U DIM2 Double)
stepPosition h qs ms ps = do
  do newQs <- computeP $ qs +^ (ps *^ h2 /^ ms2)
     return newQs
    where
      (Z :. i :. j) = extent ps

      h2  = extend (Any :. i :. j) $ fromListUnboxed Z [h]
      ms2 = extend (Any :. j) ms

stepMomentum :: forall a b c m . ( Monad m
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
stepMomentum gConst h qs ms ps =
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
    
    Z :.i :. _j :. k = extent fss
    dt2              = extend (Any :. i :. k) $ fromListUnboxed Z [h]

stepOnce :: ( Monad m
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
stepOnce gConst h ms qs ps = do
  newPs <- stepMomentum gConst h qs ms ps
  newQs <- stepPosition h qs ms newPs
  return (newQs, newPs)

repDim2to3Outer :: Source a Double =>
                   Array a DIM2 Double ->
                   Array D DIM3 Double
repDim2to3Outer a = extend (Any :. I.spaceDim) a

type Momenta   = Array U DIM2 Double
type Positions = Array U DIM2 Double
type Masses    = Array U DIM1 Double

zeroDiags :: Source a Double =>
             Array a DIM2 Double ->
             Array D DIM2 Double
zeroDiags x = traverse x id f
  where
    f _ (Z :. i :. j) | i == j    = 0.0
                      | otherwise = x!(Z :. i :. j)
                                    
zeroDiags' :: Source a Double =>
              Array a DIM3 Double ->
              Array D DIM3 Double
zeroDiags' x = traverse x id f
  where
    f _ (Z :. i :. j :. k) | i == j    = 0.0
                           | otherwise = x!(Z :. i :. j :. k)

potentialEnergy :: Double -> Masses -> Positions -> IO Double
potentialEnergy gConst ms qs = do
  ds2 <- sumP $ Repa.map (^2) $ pointDiffs qs
  let ds   = Repa.map sqrt ds2
      is   = prodPairsMasses ms
      pess = zeroDiags $ Repa.map (* (negate gConst)) $ is /^ ds
  pes <- sumP pess
  pe  <- sumP pes
  return $ head $ toList $ Repa.map (* 0.5) pe

kineticEnergy :: Masses -> Momenta -> IO Double
kineticEnergy ms ps = do
  preKes <- sumP $ ps *^ ps
  ke     <- sumP $ preKes /^ ms
  return $ head $ toList $ Repa.map (* 0.5) ke

hamiltonian :: Double -> Masses -> Momenta -> Positions -> IO Double
hamiltonian gConst ms qs ps = do
  ke <- kineticEnergy ms ps
  pe <- potentialEnergy gConst ms qs
  return $ ke + pe

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

initPs :: Array U DIM2 Double
initPs = fromListUnboxed (Z :. I.nBodiesTwoPlanets :. I.spaceDim) $ concat I.initPs

initQs :: Array U DIM2 Distance
initQs = fromListUnboxed (Z :. I.nBodiesTwoPlanets :. I.spaceDim) $ concat $ I.initQs

masses :: Array U DIM1 Mass
masses = fromListUnboxed (Z :. I.nBodiesTwoPlanets) I.massesTwoPlanets

stepN :: forall m . Monad m =>
         Int -> Double -> Double -> Masses -> Positions -> Momenta ->
         m (Positions, Momenta)
stepN n gConst dt masses = curry updaterMulti
  where
    updaterMulti = foldr (>=>) return updaters
    updaters :: [(Positions, Momenta) -> m (Positions, Momenta)]
    updaters = replicate n (uncurry (stepOnce gConst dt masses))

main :: IO ()
main = do
  hPre <- hamiltonian I.gConst masses initQs initPs
  putStrLn $ show hPre
  (newQs, newPs) <- stepN I.nStepsTwoPlanets I.gConst I.stepTwoPlanets masses initQs initPs
  hPost <- hamiltonian I.gConst masses newQs newPs
  putStrLn $ show hPost
