module Main where

import Data.List as L
import Data.Massiv.Array as A

nx :: Int
ny :: Int
nx = 100
ny = 80

eps0 :: Double
eps0 = 1.0

iterations :: Int
iterations = 10000

epsilonSim :: Array U Ix2 Double
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
{-# INLINE gradx #-}

grady :: Border Double -> Array U Ix2 Double -> Array U Ix2 Double
grady border arr = computeAs U $ mapStencil border gradienty arr
{-# INLINE grady #-}

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
    f (i :. _) x
        | i == 0 = vbottom
        | i == h - 1 = vtop
        | otherwise = x

setHBorder :: Double -> Double -> Array U Ix2 Double -> Array U Ix2 Double
setHBorder left right arr = computeAs U $ imap f arr
  where
    (_ :. w) = unSz (size arr)
    f (_ :. j) x
        | j == 0 = left
        | j == w - 1 = right
        | otherwise = x

runSim :: Double -> Int -> Int -> (Array U Ix2 Double -> Array U Ix2 Double) -> IO (Array U Ix2 Double)
runSim threshold maxIteration evalStep runfunc = runsim' 0 phi
  where
    phi = A.replicate Seq (Sz (ny :. nx)) 0.0
    runsim' n p
        | n == maxIteration = return p
        | n `mod` evalStep == 0 = do
            let p' = runfunc p
            let simerror = calcError p p'
            putStrLn $ "Iteration " ++ show n ++ " with error " ++ show simerror
            if simerror < threshold then return p' else runsim' (n + 1) p'
        | otherwise = runsim' (n + 1) (runfunc p)

calcError :: Array U Ix2 Double -> Array U Ix2 Double -> Double
calcError a1 = A.maximum' . fmap abs . A.zipWith (-) a1

main :: IO ()
main = do
    mcharge <- newMArray (Sz (ny :. nx)) 0.0 :: IO (MArray RealWorld P Ix2 Double)
    write_ mcharge (40 :. 50) (-2.0)
    write_ mcharge (20 :. 20) (-1.0)
    charge <- freeze Seq mcharge
    let charge2 = computeAs U charge

    resultphi <- runSim 1e-4 iterations 500 $ setHBorder (-0.5) 0.5 . setVBorder (-1.0) 1.0 . step charge2 epsilonSim

    let ex = gradx Wrap resultphi
        ey = grady Wrap resultphi

    write2csv "phi.csv" resultphi
    write2csv "charge.csv" $ computeAs U charge
    write2csv "epsilon.csv" $ computeAs U epsilonSim
    write2csv "Ex.csv" ex
    write2csv "Ey.csv" ey
