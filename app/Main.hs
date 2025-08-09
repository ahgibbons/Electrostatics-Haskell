module Main where

import Data.List as L
import Data.Massiv.Array as A

nx :: Int
ny :: Int
nx = 100
ny = 80

iterations :: Int
iterations = 10000

epsilonSim :: Array U Ix2 Double
-- epsilonSim = makeArrayR U Seq (Sz (ny :. nx)) (\(i :. _) -> if (i > 30 && i < 50) then 2.0 else 1.0)
epsilonSim = A.replicate Seq (Sz (ny :. nx)) 1.0

laplacianStencil :: (Fractional a) => Stencil Ix2 a a
laplacianStencil = makeStencil (Sz (3 :. 3)) (1 :. 1) $ \get ->
    0.25
        * ( get ((-1) :. 0)
                + get (1 :. 0)
                + get (0 :. 1)
                + get (0 :. (-1))
          )

gradienty :: (Fractional a) => Stencil Ix2 a a
gradienty = makeStencil (Sz (3 :. 3)) (1 :. 1) $ \get ->
    0.5 * (get (1 :. 0) - get ((-1) :. 0))
{-# INLINE gradientx #-}
gradientx :: (Fractional a) => Stencil Ix2 a a
gradientx = makeStencil (Sz (3 :. 3)) (1 :. 1) $ \get ->
    0.5 * (get (0 :. 1) - get (0 :. -1))
{-# INLINE gradienty #-}

gradx :: Border Double -> Array U Ix2 Double -> Array U Ix2 Double
gradx border arr = computeAs U $ mapStencil border gradientx arr

grady :: Border Double -> Array U Ix2 Double -> Array U Ix2 Double
grady border arr = computeAs U $ mapStencil border gradienty arr

eps0 :: Double
eps0 = 1.0

write2csv :: FilePath -> Array U Ix2 Double -> IO ()
write2csv fp arr = do
    let dtext =
            intercalate "\n"
                . L.map (intercalate "," . L.map show)
                . A.toLists
                $ arr
    writeFile fp dtext

step ::
    Array U Ix2 Double ->
    Array U Ix2 Double ->
    Array U Ix2 Double ->
    Array U Ix2 Double
step charge epsilon phi =
    computeAs U $
        A.zipWith3 (\a b c -> a + b + c) charge_term lap_term divterm
  where
    f = 0.25 / eps0

    charge_term = A.map (* f) $ A.zipWith (/) charge epsilon
    lap_term = computeAs U $ mapStencil Wrap laplacianStencil phi
    divterm =
        A.zipWith
            (+)
            (A.zipWith (*) (gradx Wrap phi) (gradx Wrap epsilon))
            (A.zipWith (*) (grady Wrap phi) (grady Wrap epsilon))

setVBorder :: Double -> Double -> Array U Ix2 Double -> Array U Ix2 Double
setVBorder vtop vbottom arr = computeAs U $ imap f arr
  where
    (h :. _) = unSz (size arr)
    f (i :. _) x =
        if i == 0
            then vbottom
            else if i == h - 1 then vtop else x
setHBorder :: Double -> Double -> Array U Ix2 Double -> Array U Ix2 Double
setHBorder left right arr = computeAs U $ imap f arr
  where
    (_ :. w) = unSz (size arr)
    f (_ :. j) x =
        if j == 0
            then left
            else if j == w - 1 then right else x

main :: IO ()
main = do
    mcharge <-
        makeMArray
            Seq
            (Sz2 ny nx)
            ( \(i :. j) ->
                return (if i == 10 && j == 30 then 2.0 else 0.0)
            ) ::
            IO (MArray RealWorld P Ix2 Double)
    write_ mcharge (40 :. 50) (-2.0)
    charge <- freeze Seq mcharge
    let charge2 = computeAs U charge
    write2csv "charge.csv" $ computeAs U charge
    write2csv "epsilon.csv" $ computeAs U epsilonSim
    let phi0' = makeArray Seq (Sz (ny :. nx)) (const 0.0) :: Array D Ix2 Double
    let phi50 = iterate (setHBorder (-1.0) (1.0) . step charge2 epsilonSim) (computeAs U phi0') !! iterations
    write2csv "phi.csv" phi50

    let ex = gradx Wrap phi50
        ey = grady Wrap phi50
    write2csv "Ex.csv" ex
    write2csv "Ey.csv" ey

    putStrLn $ "Finished running " ++ show iterations ++ " iterations"
