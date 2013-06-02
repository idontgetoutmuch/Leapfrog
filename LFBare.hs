{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE FlexibleContexts          #-}
{-# LANGUAGE TypeOperators             #-}
{-# LANGUAGE ScopedTypeVariables       #-}
{-# LANGUAGE RankNTypes                #-}

import Data.Array.Repa as Repa hiding ((++), zipWith)


type Velocities = Array U DIM2 Double
type Positions  = Array U DIM2 Double
type Masses     = Array U DIM1 Double

spaceDim = 3

gConst = 6.67384e-11  -- gravitational constant
nBodies = 3           -- number of bodies
eps2 = 0.25           -- softening constant squared
dt = 864000.0         -- timestep

sunMass     = 1.9889e30
jupiterMass = 1.8986e27
earthMass   = 5.9722e24

jupiterPerihelion = 7.405736e11
earthPerihelion = 1.470983e11

jupiterV = (-1.0965244901087316e02, -1.3710001990210707e04, 0.0)
jupiterR = (negate jupiterPerihelion, 0.0, 0.0)

earthV = (2.694354528161541e03, 3.016946927465788e04, 0.0)
earthR = (earthPerihelion, 0.0, 0.0)

sunV = (0.0, 0.0, 0.0)
sunR = (0.0, 0.0, 0.0)

initVs = fromListUnboxed (Z :. nBodies :. spaceDim) xs
  where
    xs = [ earthX,   earthY,   earthZ
         , jupiterX, jupiterY, jupiterZ
         , sunX,     sunY,     sunZ
         ]
    (earthX,   earthY,   earthZ)   = earthV
    (jupiterX, jupiterY, jupiterZ) = jupiterV
    (sunX,     sunY,     sunZ)     = sunV

initRs :: Positions
initRs = fromListUnboxed (Z :. nBodies :. spaceDim) xs
  where
    xs = [ earthX,   earthY,   earthZ
         , jupiterX, jupiterY, jupiterZ
         , sunX,     sunY,     sunZ
         ]
    (earthX,   earthY,   earthZ)   = earthR
    (jupiterX, jupiterY, jupiterZ) = jupiterR
    (sunX,     sunY,     sunZ)     = sunR

masses :: Masses
masses = fromListUnboxed (Z :. nBodies) xs
  where
    xs = [ earthMass
         , jupiterMass
         , sunMass
         ]

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

replicateRows :: Source a Double =>
                 Array a DIM2 Double ->
                 Array D DIM3 Double
replicateRows a = extend (Any :. i :. All) a
  where (Z :. i :. _j) = extent a

transposeOuter :: Source a Double =>
              Array a DIM3 Double ->
              Array D DIM3 Double
transposeOuter ps = backpermute (f e) f ps
  where
    e = extent ps
    f (Z :. i :. i' :. j) = Z :. i' :. i :. j

pointDiffs :: Source a Double =>
              Array a DIM2 Double ->
              Array D DIM3 Double
pointDiffs ps = qs -^ (transposeOuter qs)
  where qs = replicateRows ps

repDim2to3Outer :: Source a Double =>
                   Array a DIM2 Double ->
                   Array D DIM3 Double
repDim2to3Outer a = undefined -- extend (Any :. spaceDim) a

forces :: (Source a Double, Source b Double) =>
          Array a DIM2 Double ->
          Array b DIM1 Double ->
          Array D DIM3 Double
forces qs ms = fs *^ is
  where ds = repDim2to3Outer $
             Repa.map ((^3) . sqrt . (+ eps2)) $
             foldS (+) 0 $
             Repa.map (^2) $
             ps
        is = repDim2to3Outer $ prodPairsMasses ms
        fs = Repa.map (* (negate gConst)) $
             ps /^ ds
        ps = pointDiffs qs

stepVelocity :: ( Source a Double
                , Source b Double
                , Source c Double
                )  =>
                Array a DIM2 Double ->
                Array b DIM2 Double ->
                Array c DIM1 Double ->
                Array D DIM2 Double
stepVelocity vs fs ms = vs +^ (fs *^ dt2 /^ ms2)
  where
    ms2 :: Array D DIM2 Double
    ms2 = extend (Any :. j) ms
    dt2 :: Array D DIM2 Double
    dt2 = extend (Any :. i :. j) $ fromListUnboxed Z [dt]
    (Z :. i :. j) = extent fs

stepPosition :: ( Source a Double
                , Source b Double
                ) =>
                Array a DIM2 Double ->
                Array b DIM2 Double ->
                Array D DIM2 Double
stepPosition xs vs = xs +^ (vs *^ dt2)
  where
    dt2 :: Array D DIM2 Double
    dt2 = extend (Any :. i :. j) $ fromListUnboxed Z [dt]
    (Z :. i :. j) = extent vs

stepOnce :: Monad m =>
            Masses -> Positions -> Velocities ->
            m (Positions, Velocities)
stepOnce masses rs vs = do
  fs    <- sumP $ transpose $ forces rs masses
  newVs <- computeP $ stepVelocity vs fs masses
  newRs <- computeP $ stepPosition rs newVs
  return (newRs, newVs)

stepN :: forall m . Monad m =>
         Int -> Masses -> Positions -> Velocities ->
         m (Positions, Velocities)
stepN 0 masses rs vs = return (rs, vs)
stepN n masses rs vs = do
   (rs', vs') <- stepOnce masses rs vs
   stepN (n - 1) masses rs' vs'

nSteps :: Int
nSteps = 10000000

main :: IO ()
main = do
  rsVs <- stepN nSteps masses initRs initVs
  putStrLn $ show rsVs
