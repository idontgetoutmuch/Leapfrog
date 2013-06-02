import Numeric.Units.Dimensional.Prelude
import qualified Prelude
import qualified Data.Vector as V

import Text.PrettyPrint
import Debug.Trace

gConst   = 6.674e-11 *~ (newton * (meter / kilo gram) * (meter / kilo gram))

mSun     = 1.9889e30 *~ kilo gram
mJupiter = 1.8986e27 *~ kilo gram
mEarth   = 5.9722e24 *~ kilo gram

dEarth   = 1.4960e6 *~ kilo meter
dSun     = 0.0000e6 *~ kilo meter
dJupiter = 7.7855e6 *~ kilo meter

dt = 10 *~ day

forcesV :: Int -> V.Vector (V.Vector (Quantity DLength Double)) -> V.Vector (Quantity DMass Double) -> V.Vector (V.Vector (Unit DForce Double))
forcesV ix positions masses = undefined
  where
    -- Get the i-th row
    v :: V.Vector (Quantity DLength Double)
    v = positions V.! ix

    m :: Quantity DMass Double
    m = masses V.! ix

    -- Calculate the pointwise distances between the i-th particle
    -- and the rest of the particles
    ds :: V.Vector (V.Vector (Quantity DLength Double))
    ds = V.zipWith (V.zipWith (-)) positions (V.replicate (V.length positions) v)

    -- Calculate the sum of the squares of the distances between the
    -- i-th particle and the rest of the particles and "normalize"
    dsqs :: V.Vector (Quantity DLength Double)
    dsqs = V.map f ds
           where
             f :: V.Vector (Quantity DLength Double) -> (Quantity DLength Double)
             f x = undefined -- {- (**(3/2)) -} $ undefined {- V.sum -} $ V.map (^2) x

    -- -- Calculate the forces on the i-th particle
    -- fs :: V.Vector (V.Vector Double)
    -- fs = V.generate (V.length ds) f
    --      where
    --        f :: Int -> V.Vector Double
    --        f i | i == ix   = V.fromList $ take 3 $ repeat 0.0
    --            | otherwise = V.map g (ds V.! i)
    --            where
    --              g x = gConst * x * (m * (ms V.! i) / (dsqs V.! i))
