
Forget about Newton and suppose you are told that the way to do
mechanics is to write down the total energy of the system in which you
are interested and then apply Hamilton's equations.

Consider a mass of $m$ attached to a light rod of length $l$ which is
attached to a point from which it can swing freely in a plane. Then
the kinetic energy is:

$$
\frac{1}{2}mv^2 = \frac{1}{2}ml^2\dot{\theta}^2
$$

and the potential energy (taking this to be 0 at $\theta = 0$) is:

$$
mgl(1 - \cos\theta)
$$

Thus the Hamiltonian is:

$$
\cal{H} = \frac{1}{2}ml^2\dot{\theta}^2 + mgl(1 - \cos\theta)
$$

Let us set the generalized momentum

$$
p = \frac{\partial\cal{L}}{\partial\dot{\theta}} = ml^2\dot{\theta}
$$

Then we can re-write the Hamiltonian as:

$$
\cal{H} = \frac{p^2}{2ml^2} + mgl(1 - \cos\theta)
$$

Applying Hamilton's equations we obtain

$$
\begin{aligned}
\dot{\theta} &=  \frac{\partial\cal{H}}{\partial p}      = \frac{p}{ml^2} \\
\dot{p}      &= -\frac{\partial\cal{H}}{\partial \theta} = -mgl\sin\theta
\end{aligned}
$$

Differentiating the first equation with respect to time we then obtain
the familiar equation describing the motion of a simple pendulum.

$$
\ddot{\theta} = \frac{\dot{p}}{ml^2} = \frac{-mgl\sin\theta}{ml^2} = -\frac{g}{l}\sin\theta
$$

Now we would like to calculate the pendulum's position and velocity at
a given time. The obvious starting place is to use the explicit Euler
method.

$$
\begin{aligned}
\theta_{n+1} &=  h\frac{p_n}{ml^2} \\
p_{n+1}      &= -hmgl\sin\theta_n
\end{aligned}
$$

> module Symplectic (dia', main) where

> import Diagrams.Backend.Cairo.CmdLine
> import Diagrams.Prelude
>
> import Text.Printf

> stepMomentumEE :: Double -> Double -> Double -> Double -> Double
> stepMomentumEE m l p q = p -  h * m * g * l * sin q

> stepPositionEE :: Double -> Double -> Double -> Double -> Double
> stepPositionEE m l p q = q + h * p / (m * l^2)

> stepOnceEE :: Double -> Double -> Double -> Double -> (Double, Double)
> stepOnceEE m l p q = (newP, newQ)
>   where
>     newP = stepMomentumEE m l p q
>     newQ = stepPositionEE m l p q

> runEE :: Double -> Double -> [(Double, Double)]
> runEE initP initTheta = iterate (uncurry (stepOnceEE m l)) (initP, initTheta)

The symplectic Euler method:

$$
\begin{aligned}
p_{n+1} &= p_n - h\nabla_q H(p_{n+1}, q_n) \\
q_{n+1} &= q_n + h\nabla_p H(p_{n+1}, q_n)
\end{aligned}
$$

We check that this really is symplectic. First suppose we have two functions:

$$
\begin{aligned}
x &= u - f(x,v) \\
y &= v + g(x,v) \\
\end{aligned}
$$

Then we can find partial derivatives:

$$
\begin{aligned}
dx &= du - \frac{\partial f}{\partial x}dx - \frac{\partial f}{\partial v}dv \\
dy &= dv + \frac{\partial g}{\partial x}dx + \frac{\partial g}{\partial v} dv \\
\end{aligned}
$$

$$
\begin{aligned}
\frac{\partial x}{\partial u} &= 1 - \frac{\partial f}{\partial x}\frac{\partial x}{\partial u} \\
\frac{\partial x}{\partial v} &= -\frac{\partial f}{\partial x}\frac{\partial x}{\partial v} -\frac{\partial f}{\partial v} \\
\frac{\partial y}{\partial u} &= \frac{\partial g}{\partial x}\frac{\partial x}{\partial u} \\
\frac{\partial y}{\partial v} &= 1 + \frac{\partial g}{\partial x}\frac{\partial x}{\partial v} + \frac{\partial g}{\partial v}
\end{aligned}
$$

Re-arranging:

$$
\begin{aligned}
\frac{\partial x}{\partial u}(1 + \frac{\partial f}{\partial x}) &= 1 \\
\frac{\partial x}{\partial v}(1 + \frac{\partial f}{\partial x}) &= -\frac{\partial f}{\partial v} \\
\frac{\partial y}{\partial u} -\frac{\partial g}{\partial x}\frac{\partial x}{\partial u} &= 0 \\
\frac{\partial y}{\partial v} - \frac{\partial g}{\partial x}\frac{\partial x}{\partial v} &= 1 + \frac{\partial g}{\partial v}
\end{aligned}
$$

Pulling everything together in matrix form:

$$
\begin{bmatrix}
1 + \frac{\partial f}{\partial x} & 0 \\
-\frac{\partial g}{\partial x} & 1
\end{bmatrix}
\,
\begin{bmatrix}
\frac{\partial x}{\partial u} & \frac{\partial x}{\partial v} \\
\frac{\partial y}{\partial u} & \frac{\partial y}{\partial v}
\end{bmatrix}
=
\begin{bmatrix}
1 & -\frac{\partial f}{\partial v} \\
0 & 1 + \frac{\partial g}{\partial v}
\end{bmatrix}
$$


Now let's implement the symplectic Euler method for it:

$$
\begin{aligned}
p_{n+1} = p_n - hmgl\sin\theta_n \\
\theta_{n+1} = \theta_n + \frac{hp_{n+1}}{2ml^2}
\end{aligned}
$$

{-# OPTIONS_GHC -Wall                     #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing  #-}
{-# OPTIONS_GHC -fno-warn-type-defaults   #-}

> {-# LANGUAGE TupleSections                #-}
> {-# LANGUAGE NoMonomorphismRestriction    #-}


> h, m, l, g :: Double
> h = 0.01  -- Seconds
> l = 1.0  -- Metres
> m = 1.0  -- Kilograms
> g = 9.81 -- Metres * Seconds^-2

> stepMomentum :: Double -> Double -> Double -> Double -> Double
> stepMomentum m l p q = p -  h * m * g * l * sin q

> stepPosition :: Double -> Double -> Double -> Double -> Double
> stepPosition m l p q = q + h * p / (m * l^2)

> stepOnce :: Double -> Double -> Double -> Double -> (Double, Double)
> stepOnce m l p q = (newP, newQ)
>   where
>     newP = stepMomentum m l p q
>     newQ = stepPosition m l newP q

> initTheta, initThetaDot, initP :: Double
> initTheta    = 0.0
> initThetaDot = 1.7
> initP        = m * l^2 * initThetaDot

> runSE :: Double -> Double -> [(Double, Double)]
> runSE initP initTheta = iterate (uncurry (stepOnce m l)) (initP, initTheta)

> type DiagramC = Diagram Cairo R2

> tickSize :: Double
> tickSize   = 0.1
> cellColour0 = red  `withOpacity` 0.5
> cellColourEE0 = magenta `withOpacity` 0.5
> cellColour1 = blue `withOpacity` 0.5
> cellColour2 = green  `withOpacity` 0.5
> cellColour3 = yellow `withOpacity` 0.5


> gridLineWidth :: Double
> gridLineWidth = 0.001
>
> fSize :: Double
> fSize = 0.02
>
> cSize :: Double
> cSize = 0.01


> background :: DiagramC
> background = rect 1.1 1.1 # translate (r2 (0.5, 0.5))

> test :: Double -> [(AlphaColour Double, [(Double, Double)])] -> DiagramC
> test tickSize uss =
>   ticks  [0.0, tickSize..1.0] <>
>   ticksY [0.0, tickSize..1.0] <>
>   mconcat (zipWith hist clrs (normalise' zss)) <>
>   background
>   where
>     zss  = map snd uss
>     clrs = map fst uss

> hist :: AlphaColour Double -> [(Double, Double)] -> DiagramC
> hist cellColour xs = position $ hist' where
>     hist' = zip (map p2 xs) (repeat endpt)
>     tSize = 0.01
>     endpt = circle tSize # fcA cellColour  # lw 0

> normalise :: [(Double, Double)] -> [(Double, Double)]
> normalise zs = zip xs' ys'
>   where
>     xs = map fst zs
>     ys = map snd zs
>     ysmax = maximum ys
>     ysmin = minimum ys
>     xsmax = maximum xs
>     xsmin = minimum xs
>     xs' = map (\x -> (x - xsmin) / (xsmax - xsmin)) xs
>     ys' = map (\y -> (y - ysmin) / (ysmax - ysmin)) ys

> normalise' :: [[(Double, Double)]] -> [[(Double, Double)]]
> normalise' zss = zipWith zip xss' yss'
>   where
>     xss = map (map fst) zss
>     yss = map (map snd) zss
>     ysmax = maximum $ map maximum yss
>     ysmin = minimum $ map minimum yss
>     xsmax = maximum $ map maximum xss
>     xsmin = minimum $ map minimum xss
>     xss' = map (map (\x -> (x - xsmin) / (xsmax - xsmin))) xss
>     yss' = map (map (\y -> (y - ysmin) / (ysmax - ysmin))) yss

> ticks :: [Double] -> DiagramC
> ticks xs = (mconcat $ map tick xs)  <> line
>   where
>     maxX   = maximum xs
>     line   = fromOffsets [r2 (maxX, 0)]
>     tSize  = maxX / 100
>     tick x = endpt # translate tickShift
>       where
>         tickShift = r2 (x, 0)
>         endpt     = topLeftText (printf "%.2f" x) # fontSize (tSize * 2) <>
>                     circle tSize # fc red  # lw 0

> ticksY :: [Double] -> DiagramC
> ticksY xs = (mconcat $ Prelude.map tick xs)  <> line
>   where
>     maxX   = maximum xs
>     line   = fromOffsets [r2 (0, maxX)] # lw gridLineWidth
>     tick x = endpt # translate tickShift
>       where
>         tickShift = r2 (0, x)
>         endpt     = myText (printf "%.2f" x) # fontSize fSize <>
>                     circle cSize # fc red # lw 0
>         myText = alignedText 1.0 0.5

> bls   = runSE initP         initTheta
> blsEE = runEE initP         initTheta
> brs   = runSE (initP + 1.0) initTheta
> trs   = runSE (initP + 1.0) (initTheta + 1.0)
> tls   = runSE initP         (initTheta + 1.0)
>
> areaParGram (x1, y1) (x2, y2) = (x2 - x1) * (y2 - y1)
> areas = zipWith areaParGram bls trs
>
> hamiltonian :: Double -> Double -> Double -> Double -> Double
> hamiltonian m l p q = (p^2 / (2 * m * l^2)) + (m * g * l * (1 - cos q))

> nPlotPoints :: Int
> nPlotPoints = 400

> dia' :: DiagramC
> dia' = test tickSize [ (cellColour0, take nPlotPoints $ bls)
>                      , (cellColourEE0, take nPlotPoints $ blsEE)
>                      , (cellColour1, take nPlotPoints $ brs)
>                      , (cellColour2, take nPlotPoints $ trs)
>                      , (cellColour3, take nPlotPoints $ tls)
>                      ]

> main :: IO ()
> main = defaultMain dia'

```{.dia width='800'}
import Symplectic
dia = dia'
```
