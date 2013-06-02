-- import Numeric.Units.Dimensional.Prelude
-- import qualified Prelude

perihelion = 7.40573600e11 -- *~ meter
aphelion   = 8.16520800e11 -- *~ meter
eccentrity = 0.048775      -- unitless
t          = 4332.59 * 60 * 60 * 24 -- *~ day
m          = 1.8986e27     -- *~ kilo gram
ms         =  1.9889e30    -- *~ kilo gram
gConst   = 6.674e-11       -- *~ (newton * (meter / kilo gram) * (meter / kilo gram))


aphelion' = perihelion * (1 + eccentrity) / (1 - eccentrity)
a         = (aphelion + perihelion) / 2

n = sqrt $ gConst * ms / a^3
t' = 2 * pi / n

t'' = sqrt $ 4 * pi^2 * a^3 / (gConst * ms)

-- r^2 * thetȧ =  n * a^2 * sqrt (1 − e^2)

thetaDotP = n * a^2 * sqrt (1 - eccentrity^2) / perihelion^2 -- radians per second

