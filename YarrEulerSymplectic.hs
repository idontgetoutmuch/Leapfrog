{-# OPTIONS_GHC -Wall                     #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing  #-}
{-# OPTIONS_GHC -fno-warn-type-defaults   #-}

{-# LANGUAGE NoMonomorphismRestriction    #-}
{-# LANGUAGE FlexibleContexts             #-}
{-# LANGUAGE ScopedTypeVariables          #-}

module Main(main) where

import qualified Initial as I

import           Data.Yarr
import           Data.Yarr.Shape
import qualified Data.Yarr.Utils.FixedVector as V
import qualified Data.Yarr.IO.List as YIO

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

initPs :: IO Momenta
initPs = YIO.fromList I.nBodiesTwoPlanets $ map V.fromList I.initPs

initQs :: IO Distances
initQs = YIO.fromList I.nBodiesTwoPlanets $ map V.fromList I.initQs

masses :: IO Masses
masses = YIO.fromList I.nBodiesTwoPlanets I.massesTwoPlanets

vZero :: VecList N3 Double
vZero = V.replicate 0

stepPosition :: Double -> Positions -> Masses -> Momenta -> IO ()
stepPosition h qs ms vs = loadS fill (dzip3 upd qs ms vs) qs
  where
    upd q m p = V.zipWith (+) q (V.map (* (h / m)) p)

stepMomentum :: Double ->
                 Double ->
                 Positions ->
                 Masses ->
                 Momenta ->
                 IO ()
stepMomentum gConst h qs ms ps = do
  fs :: Forces <- new I.nBodiesTwoPlanets
  let forceBetween i pos1 mass1 j
        | i == j = return vZero
        | otherwise = do
          pos2 <- qs `index` j
          mass2 <- ms `index` j
          let deltas = V.zipWith (-) pos1 pos2
              dist2  = V.sum $ V.map (^2) deltas
              a = 1.0 / dist2
              b = (negate gConst) * mass1 * mass2 * a * (sqrt a)
          return $ V.map (* b) deltas
      forceAdd :: Int -> Int -> Force -> IO ()
      forceAdd i _ f = do
        f0 <- fs `index` i
        write fs i (V.zipWith (+) f0 f)
      force i pos = do
        mass <- ms `index` i
        fill (forceBetween i pos mass) (forceAdd i) 0 I.nBodiesTwoPlanets
      upd momentum force =
        V.zipWith (+) momentum (V.map (\f -> f * h) force)
  fill (index qs) force 0 I.nBodiesTwoPlanets
  loadS fill (dzip2 upd ps fs) ps

stepOnce :: Double -> Double -> Masses -> Positions -> Momenta -> IO ()
stepOnce gConst h ms qs ps = do
  stepMomentum gConst h qs ms ps
  stepPosition h qs ms ps

hamiltonianP :: Double -> Masses -> Momenta -> Positions -> IO Double
hamiltonianP gConst ms qs ps = do
  let preKes = dmap V.sum $ dzip2 (V.zipWith (*)) ps ps
      ke     = dzip2 (/) preKes ms
  undefined
--   preKes <- sumP $ ps *^ ps
--   ke     <- sumP $ preKes /^ ms
--
--   ds2 <- sumP $ Repa.map (^2) $ pointDiffs qs
--   let ds   = Repa.map sqrt ds2
--       is   = prodPairsMasses ms
--       pess = zeroDiags $ Repa.map (* (negate gConst)) $ is /^ ds
--   pes <- sumP pess
--   pe  <- sumP pes
--   te :: Array U DIM0 Double <- computeP $ ke +^ pe
--   return $ head $ toList $ Repa.map (* 0.5) te

main :: IO ()
main = do
  ms <- masses
  ps <- initPs
  qs <- initQs
  fill (\_ -> return ())
       (\_ _ -> stepOnce I.gConst I.stepTwoPlanets ms qs ps)
       0 I.nStepsTwoPlanets
  putStrLn "New qs ps yarr"
  psList <- YIO.toList ps
  putStrLn $ show psList
  qsList <- YIO.toList qs
  putStrLn $ show qsList
