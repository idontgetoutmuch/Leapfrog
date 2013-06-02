import Data.Yarr
import Data.Yarr.Repr.Delayed
import Data.Yarr.Shape
import Data.Yarr.Utils.FixedVector as V
import Data.Yarr.IO.List as Y

type Distance = VecList N3 Double
type Position = VecList N3 Double
type Force = VecList N3 Double
type Speed = VecList N3 Double
type Mass = Double

type Array = UArray F L Dim1

type Velocities = Array Speed
type Positions = Array Position
type Distances = Array Distance
type Forces = Array Force
type Masses = Array Mass

vZero :: VecList N3 Double
vZero = V.replicate 0

gConst :: Double
gConst = 6.67384e-11        -- gravitational constant

nBodies :: Int
nBodies = 3                -- number of bodies

eps2, dt :: Double
eps2 = 0.25                 -- softening constant squared
                            -- FIXME: Probably far too small
dt = 864000.0               -- timestep

sunMass, jupiterMass, earthMass :: Mass
sunMass     = 1.9889e30
jupiterMass = 1.8986e27
earthMass   = 5.9722e24

jupiterPerihelion :: Double
jupiterPerihelion = 7.405736e11

earthPerihelion :: Double
earthPerihelion = 1.470983e11

jupiterV :: Speed
jupiterV = vl_3 (-1.0965244901087316e02) (-1.3710001990210707e04) 0

jupiterR :: Distance
jupiterR = vl_3 (negate jupiterPerihelion) 0 0

earthV :: Speed
earthV = vl_3 2.694354528161541e03 3.016946927465788e04 0

earthR :: Distance
earthR = vl_3 earthPerihelion 0 0

sunV :: Speed
sunV = vZero

sunR :: Distance
sunR = vZero

initVs :: IO Velocities
initVs = Y.fromList nBodies [earthV, jupiterV, sunV]

initRs :: IO Distances
initRs = Y.fromList nBodies [earthR, jupiterR, sunR]

initMasses :: IO Masses
initMasses = Y.fromList nBodies [earthMass, jupiterMass, sunMass]

forces :: Positions -> Masses -> Forces -> IO ()
forces ps ms fs = do
    fill (\_ -> return vZero) (write fs) 0 nBodies
    let forceBetween i pos1 mass1 j
            | i == j = return vZero
            | otherwise = do
                pos2 <- ps `index` j
                mass2 <- ms `index` j
                let deltas = V.zipWith (-) pos1 pos2
                    dist2 = V.sum $ V.map (^ 2) deltas
                    a = 1.0 / (dist2 + eps2)
                    b = (negate gConst) * mass1 * mass2 * a * (sqrt a)
                return $ V.map (* b) deltas
        forceAdd i _ f = do
            f0 <- fs `index` i
            write fs i (V.zipWith (+) f0 f)
        force i pos = do
            mass <- ms `index` i
            fill (forceBetween i pos mass) (forceAdd i) 0 nBodies
    fill (index ps) force 0 nBodies

stepVelocity :: Velocities -> Forces -> Masses -> IO ()
stepVelocity vs fs ms = loadS fill (dzip3 upd vs fs ms) vs
  where upd speed force mass =
            V.zipWith (+) speed (V.map (\f -> f * dt / mass) force)

stepPosition :: Positions -> Velocities -> IO ()
stepPosition ps vs = loadS fill (dzip2 upd ps vs) ps
  where upd pos speed = V.zipWith (+) pos (V.map (* dt) speed)

stepOnce :: Masses -> Positions -> Velocities -> Forces -> IO ()
stepOnce ms ps vs fs = do
  forces ps ms fs
  stepVelocity vs fs ms
  stepPosition ps vs

nSteps :: Int
nSteps = 1000000

main :: IO ()
main = do
    vs <- initVs
    ps <- initRs
    ms <- initMasses
    fs <- new nBodies
    fill (\_ -> return ()) (\_ _ -> stepOnce ms ps vs fs) 0 nSteps
    posList <- Y.toList ps
    speedList <- Y.toList vs
    putStrLn $ show (posList, speedList)
