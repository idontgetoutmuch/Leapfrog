{-# OPTIONS_GHC -Wall                     #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing  #-}
{-# OPTIONS_GHC -fno-warn-type-defaults   #-}

{-# LANGUAGE NoMonomorphismRestriction    #-}
{-# LANGUAGE FlexibleContexts             #-}
{-# LANGUAGE ScopedTypeVariables          #-}

module Main(main) where

import qualified Data.Yarr as Y
import           Data.Yarr (loadS, dzip2, dzip3, F, L)
import           Data.Yarr.Repr.Delayed (UArray)
import           Data.Yarr.Shape (fill, Dim1)
import qualified Data.Yarr.Utils.FixedVector as V
import           Data.Yarr.Utils.FixedVector (VecList, N3)
import qualified Data.Yarr.IO.List as YIO

import Data.Array.Repa hiding ((++), zipWith, Array)
import qualified Data.Array.Repa as Repa

import qualified Data.List as L

type Distance = VecList N3 Double
type Position = VecList N3 Double
type Force = VecList N3 Double
type Speed = VecList N3 Double
type Mass = Double

type Array = UArray F L Dim1

type Momenta = Array Speed -- A lie
type Positions = Array Position
type Distances = Array Distance
type Forces = Array Force
type Masses = Array Mass

spaceDim :: Int
spaceDim = 3

vZero :: VecList N3 Double
vZero = V.replicate 0

gConst :: Double
gConst = 6.67384e-11

dt :: Double
dt = 864000.0

nBodies :: Int
nBodies = 3

sunMass, jupiterMass, earthMass :: Mass
sunMass     = 1.9889e30
jupiterMass = 1.8986e27
earthMass   = 5.9722e24

jupiterPerihelion :: Double
jupiterPerihelion = 7.405736e11

earthPerihelion :: Double
earthPerihelion = 1.470983e11

jupiterV :: Speed
jupiterV = V.vl_3 (-1.0965244901087316e02) (-1.3710001990210707e04) 0

jupiterR :: Distance
jupiterR = V.vl_3 (negate jupiterPerihelion) 0 0

earthV :: Speed
earthV = V.vl_3 2.694354528161541e03 3.016946927465788e04 0

earthR :: Distance
earthR = V.vl_3 earthPerihelion 0 0

sunV :: Speed
sunV = vZero

sunR :: Distance
sunR = vZero

initMassesL :: [Double]
initMassesL = [earthMass, jupiterMass, sunMass]

initPs :: IO Momenta
initPs = YIO.fromList nBodies ps
  where
    ps = zipWith (\m v -> V.map (*m) v) initMassesL vs
    vs = [earthV, jupiterV, sunV]

initQs :: IO Distances
initQs = YIO.fromList nBodies [earthR, jupiterR, sunR]

initMs :: IO Masses
initMs = YIO.fromList nBodies initMassesL

stepPositionY :: Double -> Positions -> Masses -> Momenta -> IO ()
stepPositionY h qs ms vs = loadS fill (dzip3 upd qs ms vs) qs
  where
    upd q m p = V.zipWith (+) q (V.map (* (h / m)) p)

stepMomentumY :: Double ->
                 Double ->
                 Positions ->
                 Masses ->
                 Momenta ->
                 IO ()
stepMomentumY gConst h qs ms ps = do
  fs :: Forces <- Y.new nBodies
  let forceBetween i pos1 mass1 j
        | i == j = return vZero
        | otherwise = do
          pos2 <- qs `Y.index` j
          mass2 <- ms `Y.index` j
          let deltas = V.zipWith (-) pos1 pos2
              dist2  = V.sum $ V.map (^2) deltas
              a = 1.0 / dist2
              b = (negate gConst) * mass1 * mass2 * a * (sqrt a)
          return $ V.map (* b) deltas
      forceAdd :: Int -> Int -> Force -> IO ()
      forceAdd i _ f = do
        f0 <- fs `Y.index` i
        Y.write fs i (V.zipWith (+) f0 f)
      force i pos = do
        mass <- ms `Y.index` i
        fill (forceBetween i pos mass) (forceAdd i) 0 nBodies
      upd momentum force =
        V.zipWith (+) momentum (V.map (\f -> f * h) force)
  fill (Y.index qs) force 0 nBodies
  loadS fill (dzip2 upd ps fs) ps

stepPositionP :: forall a b c m . ( Monad m
                 , Source a Double
                 , Source b Double
                 , Source c Double
                 ) =>
                 Double ->
                 Repa.Array a DIM2 Double ->
                 Repa.Array b DIM1 Double ->
                 Repa.Array c DIM2 Double ->
                 m (Repa.Array U DIM2 Double)
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
                 Repa.Array a DIM2 Double ->
                 Repa.Array b DIM1 Double ->
                 Repa.Array c DIM2 Double ->
                 m (Repa.Array U DIM2 Double)
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
    Z :.i :. _j :. k = extent fss
    dt2              = extend (Any :. i :. k) $ fromListUnboxed Z [h]

zeroDiags' :: Source a Double =>
              Repa.Array a DIM3 Double ->
              Repa.Array D DIM3 Double
zeroDiags' x = traverse x id f
  where
    f _ (Z :. i :. j :. k) | i == j    = 0.0
                           | otherwise = x!(Z :. i :. j :. k)

repDim2to3Outer :: Source a b => Repa.Array a DIM2 b -> Repa.Array D DIM3 b
repDim2to3Outer a = extend (Any :. spaceDim) a

prodPairsMasses :: Source a Double =>
                   Repa.Array a DIM1 Double ->
                   Repa.Array D DIM2 Double
prodPairsMasses ms = ns *^ (transpose ns)
  where
    ns = repDim1To2Outer ms

transposeOuter :: Source a Double =>
              Repa.Array a DIM3 Double ->
              Repa.Array D DIM3 Double
transposeOuter qs = backpermute (f e) f qs
  where
    e = extent qs
    f (Z :. i :. i' :. j) = Z :. i' :. i :. j

pointDiffs :: Source a Double =>
              Repa.Array a DIM2 Double ->
              Repa.Array D DIM3 Double
pointDiffs qs = qss -^ (transposeOuter qss)
  where qss = replicateRows qs

repDim1To2Outer :: Source a Double =>
                   Repa.Array a DIM1 Double ->
                   Repa.Array D DIM2 Double
repDim1To2Outer a = extend (Any :. i :. All) a
  where (Z :. i) = extent a

replicateRows :: Source a Double =>
                 Repa.Array a DIM2 Double ->
                 Repa.Array D DIM3 Double
replicateRows a = extend (Any :. i :. All) a
  where (Z :. i :. _j) = extent a

main :: IO ()
main = do
  ms <- initMs
  ps <- initPs
  qs <- initQs
  psPreList <- YIO.toList ps
  qsPreList <- YIO.toList qs
  let qsRepa = fromListUnboxed (Z :. nBodies :. spaceDim) $
               concat $
               L.transpose $
               Prelude.map (\n -> Prelude.map (V.!n) qsPreList) [0..2]
  let psRepa = fromListUnboxed (Z :. nBodies :. spaceDim) $
               concat $
               L.transpose $
               Prelude.map (\n -> Prelude.map (V.!n) psPreList) [0..2]
  putStrLn "Original ps"
  putStrLn $ show psPreList
  putStrLn $ show psRepa
  putStrLn "Original qs"
  putStrLn $ show qsPreList
  putStrLn $ show qsRepa
  newP <- stepMomentumP gConst dt
          qsRepa
          (fromListUnboxed (Z :. nBodies) initMassesL)
          psRepa
  putStrLn "New ps repa"
  putStrLn $ show newP
  stepMomentumY gConst dt qs ms ps
  putStrLn "New ps yarr"
  baz <- YIO.toList ps
  putStrLn $ show baz
