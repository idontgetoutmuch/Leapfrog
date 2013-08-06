{-# OPTIONS_GHC -Wall                     #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing  #-}
{-# OPTIONS_GHC -fno-warn-type-defaults   #-}

{-# LANGUAGE NoMonomorphismRestriction    #-}
{-# LANGUAGE FlexibleContexts             #-}
{-# LANGUAGE ScopedTypeVariables          #-}

module Initial (
    spaceDim
  , gConst
  , gConstAu
  , sunMass
  , nBodiesTwoPlanets
  , stepTwoPlanets
  , nStepsTwoPlanets
  , massesTwoPlanets
  , initQs
  , initPs
  , stepOuter
  , nStepsOuter
  , nBodiesOuter
  , massesOuter
  , initQsOuter
  , initPsOuter
    -- Possibly not needed
  , earthPerihelion
  , earthEccentrity
  , earthMajRad
  , jupiterPerihelion
  , jupiterEccentrity
  , jupiterMajRad
  ) where

import Data.List.Split

gConst :: Double
gConst = 6.67384e-11

gConstAu :: Double
gConstAu = 2.95912208286e-4

nBodiesTwoPlanets :: Int
nBodiesTwoPlanets = 3

spaceDim :: Int
spaceDim = 3

nStepsTwoPlanets :: Int
nStepsTwoPlanets = 1000

stepTwoPlanets :: Double
stepTwoPlanets = 24 * 60 * 60 * 10

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

jupiterAphelion   :: Double
jupiterAphelion   = 8.165208e11

jupiterEccentrity :: Double
jupiterEccentrity = 4.8775e-2

jupiterMajRad :: Double
jupiterMajRad = (jupiterPerihelion + jupiterAphelion) / 2

earthAphelion   :: Double
earthAphelion   = 1.520982e11

earthEccentrity :: Double
earthEccentrity = 1.6711e-2

earthMajRad :: Double
earthMajRad = (earthPerihelion + earthAphelion) / 2

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

nStepsOuter :: Int
nStepsOuter = 2000000

nBodiesOuter :: Int
nBodiesOuter = length massesOuter

stepOuter :: Double
stepOuter = 100.0

massesOuter :: [Double]
massesOuter =
  [ 9.54786104043e-4
  , 2.85583733151e-4
  , 4.37273164546e-5
  , 5.17759138449e-5
  , 1.0 / 1.3e8
  , 1.00000597682
  ]

initQsOuter :: [[Double]]
initQsOuter = chunksOf 3 xs
  where
    xs = [  -3.5023653
         ,  -3.8169847
         ,  -1.5507963
         ,   9.0755314
         ,  -3.0458353
         ,  -1.6483708
         ,   8.3101420
         , -16.2901086
         ,  -7.2521278
         ,  11.4707666
         , -25.7294829
         , -10.8169456
         , -15.5387357
         , -25.2225594
         ,  -3.1902382
         ,   0.0
         ,   0.0
         ,   0.0
         ]

initPsOuter :: [[Double]]
initPsOuter = chunksOf 3 ys
  where
  ys = zipWith (*) (concat $ map (replicate 3) massesOuter) xs
  xs = [  0.00565429
       , -0.00412490
       , -0.00190589
       ,  0.00168318
       ,  0.00483525
       ,  0.00192462
       ,  0.00354178
       ,  0.00137102
       ,  0.00055029
       ,  0.00288930
       ,  0.00114527
       ,  0.00039677
       ,  0.00276725
       , -0.00170702
       , -0.00136504
       ,  0.0
       ,  0.0
       ,  0.0
       ]
