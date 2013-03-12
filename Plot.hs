{-# LANGUAGE NoMonomorphismRestriction, TupleSections #-}
module Plot
  ( plot
  , scatterplot
  , grid
  , gridSquare
  , histogram
  ) where

import Diagrams.Prelude

gridSquare = alignBL $ unitSquare # fc beige # lc white # lw 0.001

-- | Plots a function.
plot func xs = scaleX sX . scaleY sY . alignBL . position $ plotf ++ grid where
    ys = map func xs
    ysmax = fromInteger . ceiling . maximum $ ys
    ysmin = fromInteger . floor . minimum $ ys
    xsmax = fromInteger . ceiling . maximum $ xs
    xsmin = fromInteger . floor . minimum $ xs
    sX = 1 / (xsmax - xsmin + 1)
    sY = 1 / (ysmax - ysmin)
    grid = zip [p2 (x,y) | x <- [xsmin-1..xsmax-1], y <- [ysmin..ysmax-1]] (repeat gridSquare)
    plotf = zip (zipWith (\x y -> p2 (x,y)) xs ys) (repeat dot)
    dot = circle 0.1 # fc blue # lw 0.005

-- | Creates a scatteplot.
scatterplot ps = scaleX sX . scaleY sY . alignBL . position $ plotp ++ grid where
    ysmax = fromInteger . ceiling . maximum $ map snd ps
    ysmin = fromInteger . floor . minimum $ map snd ps
    xsmax = fromInteger . ceiling . maximum $ map fst ps
    xsmin = fromInteger . floor . minimum $ map fst ps
    sX = 1 / (xsmax - xsmin + 1)
    sY = 1 / (ysmax - ysmin)
    grid = zip [p2 (x,y) | x <- [xsmin-1..xsmax-1], y <- [ysmin..ysmax-1]] (repeat gridSquare)
    plotp = zip (map p2 ps) (repeat dot)
    dot = circle 0.1 # fc blue # lw 0.005

histogram ps = scaleX sX . scaleY sY . alignBL . position $ hist ++ grid where
    ysmax = fromInteger . ceiling $ maximum ps
    ysmin = fromInteger . floor $ minimum ps
    xsmax = fromIntegral $ length ps
    xsmin = 0.0
    sX = 1 / (xsmax - xsmin + 1)
    sY = 1 / (ysmax - ysmin)
    grid = zip [p2 (x,y) | x <- [xsmin-1..xsmax-1], y <- [ysmin..ysmax-1]] (repeat gridSquare)

    cell w h = alignBL $ rect w h # fc blue # lc white # lw 0.001 # opacity 0.7
    hist = zip (map p2 $ map (,0) $ map fromIntegral [0..])
               (map (cell 1) ps)

-- | Draws an empty grid.
grid m n = scaleX (1/m) . scaleY (1/n) . alignBL. position $ grid' where
    grid' = zip [p2 (x,y) | x <- [1..m], y <- [1..n]] (repeat gridSquare)