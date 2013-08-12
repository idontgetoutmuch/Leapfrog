{-# OPTIONS_GHC -Wall                     #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing  #-}
{-# OPTIONS_GHC -fno-warn-type-defaults   #-}

{-# LANGUAGE NoMonomorphismRestriction    #-}
{-# LANGUAGE FlexibleContexts             #-}
{-# LANGUAGE ScopedTypeVariables          #-}
{-# LANGUAGE GeneralizedNewtypeDeriving   #-}
{-# LANGUAGE TypeOperators                #-}
{-# LANGUAGE BangPatterns                 #-}

module RepaOnly (
    main
  ) where

import           Data.Array.Repa hiding ((++), zipWith)
import qualified Data.Array.Repa as Repa


stepOuter :: Double
stepOuter = 100.0

nStepsOuter :: Int
nStepsOuter = 2000000

spaceDim :: Int
spaceDim = 3

gConstAu :: Double
gConstAu = 2.95912208286e-4

newtype PositionP a = QP { positionP :: Array a DIM2 Double }
newtype MomentaP  a = PP { momentaP :: Array a DIM2 Double }
newtype MassP     a = MP { massP :: Array a DIM1 Double }

stepPositionP :: forall a b c m .
                 ( Monad m
                 , Source a Double
                 , Source b Double
                 , Source c Double
                 ) =>
                 Double ->
                 PositionP a ->
                 MassP b ->
                 MomentaP c ->
                 m (PositionP U)
stepPositionP !h !qs !ms !ps = do
  do newQs <- computeP $ (positionP qs) +^ ((momentaP ps) *^ h2 /^ ms2)
     return $ QP newQs
    where
      (Z :. i :. j) = extent $ momentaP ps

      h2  = extend (Any :. i :. j) $ fromListUnboxed Z [h]
      ms2 = extend (Any :. j) $ massP ms
{-# INLINE stepPositionP #-}

stepMomentumP :: forall a b c m .
                 ( Monad m
                 , Source a Double
                 , Source b Double
                 , Source c Double
                 ) =>
                 Double ->
                 Double ->
                 PositionP a ->
                 MassP b ->
                 MomentaP c ->
                 m (MomentaP U)
stepMomentumP !gConst !h !qs !ms !ps =
  do fs <- sumP $ transpose $ zeroDiags fss
     newPs <- computeP $ (momentaP ps) +^ (fs *^ dt2)
     return $ PP newPs
  where
    is = repDim2to3Outer $ prodPairsMasses $ massP ms
    qDiffs = pointDiffs $ positionP qs
    preDs = Repa.map (^3) $
            Repa.map sqrt $
            sumS $
            Repa.map (^2) $
            qDiffs
    ds    = repDim2to3Outer preDs
    preFs = Repa.map (* (negate gConst)) $
            qDiffs /^ ds
    fss = is *^ preFs
    Z :.i :. _j :. k = extent fss
    dt2              = extend (Any :. i :. k) $ fromListUnboxed Z [h]

    repDim2to3Outer !a = extend (Any :. spaceDim) a

    zeroDiags !x = traverse x id f
      where
        f _ !(Z :. i :. j :. k) | i == j    = 0.0
                                | otherwise = x!(Z :. i :. j :. k)
{-# INLINE stepMomentumP #-}

stepOnceP :: ( Monad m
             , Source a Double
             , Source b Double
             , Source c Double
             ) =>
             Double ->
             Double ->
             MassP b ->
             PositionP a ->
             MomentaP c ->
             m (PositionP U, MomentaP U)
stepOnceP !gConst !h !ms !qs !ps = do
  newPs <- stepMomentumP gConst h qs ms ps
  newQs <- stepPositionP h qs ms newPs
  return (newQs, newPs)

prodPairsMasses :: Source a Double =>
                   Array a DIM1 Double ->
                   Array D DIM2 Double
prodPairsMasses !ms = ns *^ (transpose ns)

  where
    (Z :. i) = extent ms
    ns = extend (Any :. i :. All) ms
{-# INLINE prodPairsMasses #-}

pointDiffs :: Source a Double =>
              Array a DIM2 Double ->
              Array D DIM3 Double
pointDiffs !qs = qss -^ (transposeOuter qss)
  where

    qss = replicateRows qs

    transposeOuter qs = backpermute (f e) f qs
      where
        e = extent qs
        f !(Z :. i :. i' :. j) = Z :. i' :. i :. j

    replicateRows :: Source a Double =>
                     Array a DIM2 Double ->
                     Array D DIM3 Double
    replicateRows !a = extend (Any :. i :. All) a
      where (Z :. i :. _j) = extent a
{-# INLINE pointDiffs #-}

stepN :: forall m . Monad m =>
         Int ->
         Double ->
         Double ->
         MassP U ->
         PositionP U ->
         MomentaP U ->
         m (PositionP U, MomentaP U)
stepN !n !gConst !dt !masses !qs !ps = go n qs ps
  where
    go !i !qs !ps
      | i == 0
      = return (qs, ps)
      | otherwise
      = do (qs', ps') <- stepOnceP gConst dt masses qs ps
           go (i - 1) qs' ps'
{-# INLINE stepN #-}

mosssP :: MassP U
mosssP = MP $ fromListUnboxed (Z :. (2 :: Int))
         [ 9.54786104043e-4
         , 1.00000597682
         ]

qosss :: PositionP U
qosss = QP $ fromListUnboxed ((Z :. 2) :. 3)
        [  -3.5023653,-3.8169847,-1.5507963
        ,   0.0,0.0,0.0]

posss :: MomentaP U
posss = PP $ fromListUnboxed ((Z :. 2) :. 3)
        [  5.398637520229294e-6,-3.938397200566971e-6,-1.8197172878345132e-6
        ,  0.0,0.0,0.0]

main :: IO ()
main = do
  (qsPost, psPost) <- stepN nStepsOuter gConstAu stepOuter
                      mosssP qosss posss
  putStrLn $ show $ positionP qsPost
  putStrLn $ show $ momentaP psPost
