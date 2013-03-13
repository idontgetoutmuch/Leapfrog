{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE TupleSections             #-}
{-# LANGUAGE TypeFamilies              #-}
{-# LANGUAGE FlexibleContexts          #-}

{-# OPTIONS_GHC -Wall                    #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing #-}
{-# OPTIONS_GHC -fno-warn-type-defaults  #-}

module Plot
  ( plot
  , background
  , ticksX
  , ticksY
  ) where

import Diagrams.Prelude
import Text.Printf
import Diagrams.TwoD.Text

gridLineWidth :: Double
gridLineWidth = 0.001

fSize :: Double
fSize = 0.02

cSize :: Double
cSize = 0.001

background s = rect (s * 1.2) (s * 1.2) # showOrigin

ticksX :: (Renderable (Path R2) b, Renderable Text b) =>
          [Double] -> Diagram b R2
ticksX xs = (mconcat $ Prelude.map tick xs)  <> line
  where
    maxX   = maximum xs
    line   = fromOffsets [r2 (maxX, 0)] # lw gridLineWidth
    tick x = endpt # translate tickShift
      where
        tickShift = r2 (x, 0)
        endpt     = topLeftText (printf "%.2f" x) # fontSize fSize <>
                    circle cSize # fc green # lw 0

ticksY :: (Renderable Text b, Renderable (Path R2) b) =>
          [Double] -> Diagram b R2
ticksY xs = (mconcat $ Prelude.map tick xs)  <> line
  where
    maxX   = maximum xs
    line   = fromOffsets [r2 (0, maxX)] # lw gridLineWidth
    tick x = endpt # translate tickShift
      where
        tickShift = r2 (0, x)
        endpt     = myText (printf "%.2f" x) # fontSize fSize <>
                    circle cSize # fc green # lw 0
        myText = alignedText 1.0 0.5

plot sX sY ps = mconcat qs
  where
    qs  = scaleX sX $
          scaleY sY $
          map (\x -> dot # translate (r2 x)) ps
    dot = circle 0.01 # fc blue # lw 0.005