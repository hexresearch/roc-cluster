module Main where

import Data.Cluster.ROC
import Data.Monoid
import Graphics.Gloss.Interface.Pure.Game
import Options.Applicative

import qualified Data.Foldable as F
import qualified Data.Sequence as S
import Debug.Trace

data Vec2 = Vec2 !Float !Float
  deriving (Eq, Show)

instance ClusterSpace Vec2 where
  pointZero = Vec2 0 0
  pointAdd (Vec2 x1 y1) (Vec2 x2 y2) = Vec2 (x1 + x2) (y1 + y2)
  pointScale v (Vec2 x y) = Vec2 (realToFrac v * x) (realToFrac v * y)
  -- gaus kernell with sigma = 1
  pointKernel (Vec2 x1 y1) (Vec2 x2 y2) = realToFrac $ exp ( negate $ 0.5 * ((x2 - x1) ^ 2 + (y2 - y1) ^ 2) )
  pointDistanceSquared x y = 2 - 2 * pointKernel x y
  {-# INLINE pointZero #-}
  {-# INLINE pointAdd #-}
  {-# INLINE pointScale #-}
  {-# INLINE pointKernel #-}
  {-# INLINE pointDistanceSquared #-}

data World = World {
  worldClusters  :: ROCContext Vec2
, worldWidth     :: Float
, worldHeight    :: Float
, worldNewPoints :: S.Seq Vec2
, worldOldPoints :: S.Seq Vec2
}

runClustering :: ROCConfig -> IO ()
runClustering cfg = play display bg fps world render handleEvent update
  where
    width = 800
    height = 800
    display = InWindow "ROC demo" (width, height) (100, 100)
    bg = greyN 0.9
    fps = 60
    world = World (emptyROCContext cfg) (fromIntegral width) (fromIntegral height) mempty mempty

render :: World -> Picture
render World{..} = scaleWorld $ drawOrigin <> drawNewPoints <> drawOldPoints <> drawClusters
  where
    scaleWorld = scale (worldWidth*0.5) (worldHeight*0.5)
    drawOrigin =
         line [(-0.9, 0), (0.9, 0)]
      <> line [(0, -0.9), (0, 0.9)]
    drawNewPoints = F.foldMap (\(Vec2 x y) -> translate x y $ color red $ circleSolid 0.005) worldNewPoints
    drawOldPoints = F.foldMap (\(Vec2 x y) -> translate x y $ color blue $ circleSolid 0.005) worldOldPoints
    drawClusters = F.foldMap drawCluster $ rocPrototypes worldClusters
    drawCluster p = let
      Vec2 x y = prototypeValue p
      r = realToFrac $ 0.1 / (1 + exp (negate $ prototypeWeight p) )
      in translate x y $ circle r <> line [(-0.005, 0), (0.005, 0)] <> line [(0, -0.005), (0, 0.005)]

handleEvent :: Event -> World -> World
handleEvent (EventResize (w, h)) world = world {
    worldWidth = fromIntegral w
  , worldHeight = fromIntegral h
  }
handleEvent (EventKey (MouseButton LeftButton) Down _ (x, y)) world@World{..} = world {
    worldNewPoints = worldNewPoints S.|> Vec2 (2 * x / worldWidth) (2 * y / worldHeight)
  }
handleEvent (EventKey (MouseButton RightButton) Down _ _) world@World{..} = world {
    worldNewPoints = mempty
  , worldOldPoints = worldOldPoints <> worldNewPoints
  , worldClusters = clusterize worldNewPoints worldClusters
  }
handleEvent _ world = world

update :: Float -> World -> World
update = const id

data Options = Options {
  threshold :: Double
, maxClusters :: Int
}

toROCConfig :: Options -> ROCConfig
toROCConfig Options{..} = defaultROCConfig {
    rocThreshold = threshold
  , rocMaxClusters = maxClusters
  }

configParser :: Parser ROCConfig
configParser = fmap toROCConfig $ Options
  <$> option auto (
       long "threshold"
    <> help "Cluster with lower weight will be deleted at post-process step."
    <> value 0.1
    )
  <*> option auto (
       long "maxclusters"
    <> help "Maximum count of allowed clusters."
    <> value 10
    )

main :: IO ()
main = runClustering =<< execParser opts
  where
    opts = info (configParser <**> helper)
      ( fullDesc
     <> progDesc "Run interactive demo for ROC clustering algorithm"
     <> header "roc-cluster-demo - a demo for roc-cluster-demo" )
