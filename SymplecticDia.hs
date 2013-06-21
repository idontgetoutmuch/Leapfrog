module SymplecticDia (
    DiagramC
  , test
  , tickSize
  , nPlotPoints
  , cellColourEE0
  , cellColour0
  , cellColour1
  , cellColour2
  , cellColour3
 ) where

import Diagrams.Backend.Cairo.CmdLine
import Diagrams.Prelude

import Text.Printf


type DiagramC = Diagram Cairo R2

nPlotPoints :: Int
nPlotPoints = 400

background :: DiagramC
background = rect 1.1 1.1 # translate (r2 (0.5, 0.5))

gridLineWidth :: Double
gridLineWidth = 0.001

fSize :: Double
fSize = 0.02

cSize :: Double
cSize = 0.01

tickSize :: Double
tickSize   = 0.1
cellColour0 = red  `withOpacity` 0.5
cellColourEE0 = magenta `withOpacity` 0.5
cellColour1 = blue `withOpacity` 0.5
cellColour2 = green  `withOpacity` 0.5
cellColour3 = yellow `withOpacity` 0.5

test :: Double -> [(AlphaColour Double, [(Double, Double)])] -> DiagramC
test tickSize uss =
  ticks  [0.0, tickSize..1.0] <>
  ticksY [0.0, tickSize..1.0] <>
  mconcat (zipWith hist clrs (normalise' zss)) <>
  background
  where
    zss  = map snd uss
    clrs = map fst uss

hist :: AlphaColour Double -> [(Double, Double)] -> DiagramC
hist cellColour xs = position $ hist' where
    hist' = zip (map p2 xs) (repeat endpt)
    tSize = 0.002
    endpt = circle tSize # fcA cellColour  # lw 0

normalise :: [(Double, Double)] -> [(Double, Double)]
normalise zs = zip xs' ys'
  where
    xs = map fst zs
    ys = map snd zs
    ysmax = maximum ys
    ysmin = minimum ys
    xsmax = maximum xs
    xsmin = minimum xs
    xs' = map (\x -> (x - xsmin) / (xsmax - xsmin)) xs
    ys' = map (\y -> (y - ysmin) / (ysmax - ysmin)) ys

normalise' :: [[(Double, Double)]] -> [[(Double, Double)]]
normalise' zss = zipWith zip xss' yss'
  where
    xss = map (map fst) zss
    yss = map (map snd) zss
    ysmax = maximum $ map maximum yss
    ysmin = minimum $ map minimum yss
    xsmax = maximum $ map maximum xss
    xsmin = minimum $ map minimum xss
    xss' = map (map (\x -> (x - xsmin) / (xsmax - xsmin))) xss
    yss' = map (map (\y -> (y - ysmin) / (ysmax - ysmin))) yss

ticks :: [Double] -> DiagramC
ticks xs = (mconcat $ map tick xs)  <> line
  where
    maxX   = maximum xs
    line   = fromOffsets [r2 (maxX, 0)]
    tSize  = maxX / 100
    tick x = endpt # translate tickShift
      where
        tickShift = r2 (x, 0)
        endpt     = topLeftText (printf "%.2f" x) # fontSize (tSize * 2) <>
                    circle tSize # fc red  # lw 0

ticksY :: [Double] -> DiagramC
ticksY xs = (mconcat $ Prelude.map tick xs)  <> line
  where
    maxX   = maximum xs
    line   = fromOffsets [r2 (0, maxX)] # lw gridLineWidth
    tick x = endpt # translate tickShift
      where
        tickShift = r2 (0, x)
        endpt     = myText (printf "%.2f" x) # fontSize fSize <>
                    circle cSize # fc red # lw 0
        myText = alignedText 1.0 0.5

