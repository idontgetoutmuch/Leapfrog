{-# OPTIONS_GHC -Wall                     #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing  #-}
{-# OPTIONS_GHC -fno-warn-type-defaults   #-}

{-# LANGUAGE NoMonomorphismRestriction    #-}
{-# LANGUAGE FlexibleContexts             #-}
{-# LANGUAGE ScopedTypeVariables          #-}

module Main (main) where

import qualified Initial as I

import           Data.Yarr
import           Data.Yarr.Walk
import qualified Data.Yarr.Shape as S
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
stepPosition h qs ms vs = loadS S.fill (dzip3 upd qs ms vs) qs
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
        S.fill (forceBetween i pos mass) (forceAdd i) 0 I.nBodiesTwoPlanets
      upd momentum force =
        V.zipWith (+) momentum (V.map (\f -> f * h) force)
  S.fill (index qs) force 0 I.nBodiesTwoPlanets
  loadS S.fill (dzip2 upd ps fs) ps

stepOnce :: Double -> Double -> Masses -> Positions -> Momenta -> IO ()
stepOnce gConst h ms qs ps = do
  stepMomentum gConst h qs ms ps
  stepPosition h qs ms ps

potentialEnergy :: Double -> Masses -> Positions -> IO (Array Double)
potentialEnergy gConst ms qs = do
  pes :: Array Double <- new I.nBodiesTwoPlanets
  let peOnePairParticles :: Int ->
                            Int ->
                            IO Double
      peOnePairParticles i j
        | i == j = return 0.0
        | otherwise = do
          q1 <- qs `index` i
          m1 <- ms `index` i
          q2 <- qs `index` j
          m2 <- ms `index` j
          let qDiffs = V.zipWith (-) q1 q2
              dist2  = V.sum $ V.map (^2) qDiffs
              a      = 1.0 / dist2
              b      = (negate gConst) * m1 * m2 * (sqrt a)
          return b
      peAdd i _ pe = do
        peDelta <- pes `index` i
        write pes i (pe + peDelta)
      peFn i _ = do
        S.fill (peOnePairParticles i) (peAdd i) 0 I.nBodiesTwoPlanets
  S.fill (index qs) peFn 0 I.nBodiesTwoPlanets
  return pes

hamiltonian :: Double -> Masses -> Positions -> Momenta-> IO Double
hamiltonian gConst ms qs ps = do
  let preKes = dmap V.sum $ dzip2 (V.zipWith (*)) ps ps
      kes     = dzip2 (/) preKes (delay ms)
  ke <- walk (reduceL S.foldl (+)) (return 0) kes
  pes :: Array Double <- new I.nBodiesTwoPlanets
  let peOnePairParticles :: Int ->
                            Int ->
                            IO Double
      peOnePairParticles i j
        | i == j = return 0.0
        | otherwise = do
          q1 <- qs `index` i
          m1 <- ms `index` i
          q2 <- qs `index` j
          m2 <- ms `index` j
          let qDiffs = V.zipWith (-) q1 q2
              dist2  = V.sum $ V.map (^2) qDiffs
              a      = 1.0 / dist2
              b      = (negate gConst) * m1 * m2 * (sqrt a)
              c      =  V.sum $ V.map (*b) qDiffs
          return c
      peAdd i _ pe = do
        peDelta <- pes `index` i
        write pes i (pe + peDelta)
      peFn i _ = do
        S.fill (peOnePairParticles i) (peAdd i) 0 I.nBodiesTwoPlanets
  S.fill (index qs) peFn 0 I.nBodiesTwoPlanets
  pe <- walk (reduceL S.foldl (+)) (return 0) pes
  return $ pe + ke
  
main :: IO ()
main = do
  ms <- masses
  ps <- initPs
  qs <- initQs
  potEnPre <- potentialEnergy I.gConst ms qs
  potEnPreList <- YIO.toList potEnPre
  putStrLn $ show potEnPreList
