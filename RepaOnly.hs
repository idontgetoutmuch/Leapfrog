{-# OPTIONS_GHC -Wall                     #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing  #-}
{-# OPTIONS_GHC -fno-warn-type-defaults   #-}

{-# LANGUAGE NoMonomorphismRestriction    #-}
{-# LANGUAGE FlexibleContexts             #-}
{-# LANGUAGE ScopedTypeVariables          #-}
{-# LANGUAGE GeneralizedNewtypeDeriving   #-}
{-# LANGUAGE TypeOperators                #-}

module RepaOnly (
    outerPlanets
  , main
  ) where

import           Data.Array.Repa hiding ((++), zipWith)
import qualified Data.Array.Repa as Repa

import           Control.Monad
import           Control.Monad.Identity

import qualified Initial as I


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
stepPositionP h qs ms ps = do
  do newQs <- computeP $ (positionP qs) +^ ((momentaP ps) *^ h2 /^ ms2)
     return $ QP newQs
    where
      (Z :. i :. j) = extent $ momentaP ps

      h2  = extend (Any :. i :. j) $ fromListUnboxed Z [h]
      ms2 = extend (Any :. j) $ massP ms

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
stepMomentumP gConst h qs ms ps =
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

    repDim2to3Outer a = extend (Any :. I.spaceDim) a

    zeroDiags x = traverse x id f
      where
        f _ (Z :. i :. j :. k) | i == j    = 0.0
                               | otherwise = x!(Z :. i :. j :. k)

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
stepOnceP gConst h ms qs ps = do
  newPs <- stepMomentumP gConst h qs ms ps
  newQs <- stepPositionP h qs ms newPs
  return (newQs, newPs)


kineticEnergyP :: MassP U -> MomentaP U -> IO (Array D DIM0 Double)
kineticEnergyP ms ps = do
  preKes <- sumP $ (momentaP ps) *^ (momentaP ps)
  ke     <- sumP $ preKes /^ (massP ms)
  return $ Repa.map (* 0.5) ke

potentialEnergyP :: Double ->
                    MassP U ->
                    PositionP U ->
                    IO (Array U DIM1 Double)
potentialEnergyP gConst ms qs = do
  ds2 <- sumP $ Repa.map (^2) $ pointDiffs $ positionP qs
  let ds   = Repa.map sqrt ds2
      is   = prodPairsMasses $ massP ms
      pess = zeroDiags $ Repa.map (* (0.5 * negate gConst)) $ is /^ ds
  pes <- sumP pess
  return pes

  where

    zeroDiags x = traverse x id f
      where
        f _ (Z :. i :. j) | i == j    = 0.0
                          | otherwise = x!(Z :. i :. j)

hamiltonianP :: Double ->
                MassP U ->
                PositionP U ->
                MomentaP U ->
                IO Double
hamiltonianP gConst ms qs ps = do
  ke <- kineticEnergyP ms ps
  pes <- potentialEnergyP gConst ms qs
  pe  <- sumP pes
  te :: Array U DIM0 Double <- computeP $ ke +^ pe
  return $ head $ toList te

prodPairsMasses :: Source a Double =>
                   Array a DIM1 Double ->
                   Array D DIM2 Double
prodPairsMasses ms = ns *^ (transpose ns)

  where
    (Z :. i) = extent ms
    ns = extend (Any :. i :. All) ms

pointDiffs :: Source a Double =>
              Array a DIM2 Double ->
              Array D DIM3 Double
pointDiffs qs = qss -^ (transposeOuter qss)
  where

    qss = replicateRows qs

    transposeOuter qs = backpermute (f e) f qs
      where
        e = extent qs
        f (Z :. i :. i' :. j) = Z :. i' :. i :. j

    replicateRows :: Source a Double =>
                     Array a DIM2 Double ->
                     Array D DIM3 Double
    replicateRows a = extend (Any :. i :. All) a
      where (Z :. i :. _j) = extent a

stepN :: forall m . Monad m =>
         Int ->
         Double ->
         Double ->
         MassP U ->
         PositionP U ->
         MomentaP U ->
         m (PositionP U, MomentaP U)
stepN n gConst dt masses = curry updaterMulti
  where
    updaterMulti = foldr (>=>) return updaters
    updaters = replicate n (uncurry (stepOnceP gConst dt masses))

stepNs :: Monad m =>
          Int ->
          Double ->
          Double ->
          MassP U ->
          PositionP U ->
          MomentaP U ->
          m [(PositionP U, MomentaP U)]
stepNs n gConst dt ms rs vs = do
  rsVs <- stepAux n rs vs
  return $ (rs, vs) : rsVs
  where
    stepAux 0  _  _ = return []
    stepAux n rs vs = do
      (newRs, newVs) <- stepOnceP gConst dt ms rs vs
      rsVs <- stepAux (n-1) newRs newVs
      return $ (newRs, newVs) : rsVs

mosssP :: MassP U
mosssP = MP $ fromListUnboxed (Z :. n) I.massesOuter
  where
    n = length I.massesOuter

qosss :: PositionP U
qosss = QP $ fromListUnboxed (Z :. n :. I.spaceDim) xs
  where
    xs = concat I.initQsOuter
    n  = length xs `div` I.spaceDim

posss :: MomentaP U
posss = PP $ fromListUnboxed (Z :. n :. I.spaceDim) xs
  where
    xs = concat I.initPsOuter
    n  = length xs `div` I.spaceDim

outerPlanets :: [[(Double, Double)]]
outerPlanets = runIdentity $ do
  rsVs <- stepNs 2000 I.gConstAu 100 mosssP qosss posss
  let qs = Prelude.map fst rsVs
      xxs = Prelude.map
            (\i -> Prelude.map ((!(Z :. (i :: Int) :. (0 :: Int))) .
                                positionP) qs)
            [5,0,1,2,3,4]
      xys = Prelude.map
            (\i -> Prelude.map ((!(Z :. (i :: Int) :. (1 :: Int))) .
                                positionP) qs)
            [5,0,1,2,3,4]
  return $ zipWith zip xxs xys

main :: IO ()
main = do
  hPre <- hamiltonianP I.gConstAu mosssP qosss posss
  putStrLn $ show hPre
  (qsPost, psPost) <- stepN I.nStepsOuter I.gConstAu I.stepOuter
                      mosssP qosss posss
  hPost <- hamiltonianP I.gConstAu mosssP qsPost psPost
  putStrLn $ show hPost
