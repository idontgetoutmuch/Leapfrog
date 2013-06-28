{-# OPTIONS_GHC -Wall                     #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing  #-}
{-# OPTIONS_GHC -fno-warn-type-defaults   #-}

{-# LANGUAGE NoMonomorphismRestriction    #-}
{-# LANGUAGE FlexibleContexts             #-}
{-# LANGUAGE ScopedTypeVariables          #-}

module Initial (
    spaceDim
  , gConst
  , nBodiesTwoPlanets
  , stepTwoPlanets
  , nStepsTwoPlanets
  , massesTwoPlanets
  , initQs
  , initPs
  ) where

gConst :: Double
gConst = 6.67384e-11

nBodiesTwoPlanets :: Int
nBodiesTwoPlanets = 3

spaceDim :: Int
spaceDim = 3

nStepsTwoPlanets :: Int
nStepsTwoPlanets = 3

stepTwoPlanets :: Double
stepTwoPlanets = 864000.0

sunMass, jupiterMass, earthMass :: Double
sunMass     = 1.9889e30
jupiterMass = 1.8986e27
earthMass   = 5.9722e24

massesTwoPlanets :: [Double]
massesTwoPlanets = [earthMass, jupiterMass, sunMass]

jupiterPerihelion :: Double
jupiterPerihelion = 7.405736e11

earthPerihelion :: Double
earthPerihelion = 1.470983e11

jupiterV :: [Double]
jupiterV = [-1.0965244901087316e02, -1.3710001990210707e04, 0.0]

jupiterQ :: [Double]
jupiterQ = [negate jupiterPerihelion, 0.0, 0.0]

earthV :: [Double]
earthV = [2.694354528161541e03, 3.016946927465788e04, 0.0]

earthQ :: [Double]
earthQ = [earthPerihelion, 0.0, 0.0]

sunV :: [Double]
sunV = [0.0, 0.0, 0.0]

sunQ :: [Double]
sunQ = [0.0, 0.0, 0.0]

initPs :: [[Double]]
initPs = zipWith (\m v -> map (*m) v)
         massesTwoPlanets
         [earthV, jupiterV, sunV]

initQs :: [[Double]]
initQs = [earthQ, jupiterQ, sunQ]

