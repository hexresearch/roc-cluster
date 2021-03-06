module Data.Cluster.ROCSpec(
    spec
  ) where

import Data.Cluster.ROC
import Data.List (nub)
import Test.Hspec
import Test.HUnit

-- | Config without threshold
cfg :: ROCConfig
cfg = defaultROCConfig

-- | Config with threshold
cfg' :: ROCConfig
cfg' = defaultROCConfig {
    rocThreshold = 0.1
  }

data Vec2 = Vec2 !Double !Double
  deriving (Eq, Show)

instance ClusterSpace Vec2 where
  pointZero = Vec2 0 0
  pointAdd (Vec2 x1 y1) (Vec2 x2 y2) = Vec2 (x1 + x2) (y1 + y2)
  pointScale v (Vec2 x y) = Vec2 (v*x) (v*y)
  -- gaus kernell with sigma = 1
  pointKernel (Vec2 x1 y1) (Vec2 x2 y2) = exp ( negate $ 0.5 * ((x2 - x1) ^ 2 + (y2 - y1) ^ 2) )
  pointDistanceSquared x y = 2 - 2 * pointKernel x y
  {-# INLINE pointZero #-}
  {-# INLINE pointAdd #-}
  {-# INLINE pointScale #-}
  {-# INLINE pointKernel #-}
  {-# INLINE pointDistanceSquared #-}

spec :: Spec
spec = do
  it "Handle empty clusterization" $ do
    let cntx :: ROCContext Vec2 = clusterize [] $ emptyROCContext cfg
    assertEqual "should be empty" [] $ rocPrototypes cntx
  it "Handle single clusterization" $ do
    let cntx :: ROCContext Vec2 = clusterize [Vec2 0 0] $ emptyROCContext cfg
    assertEqual "should be empty" [newPrototype (Vec2 0 0)] $ rocPrototypes cntx
  it "Handle simple clusterization" $ do
    let cntx :: ROCContext Vec2 = clusterize [Vec2 0 0, Vec2 0.1 0] $ emptyROCContext cfg
    assertBool "should be 2 clusters" $ length (nub $ rocPrototypes cntx) == 2
  it "Handle simple clusterization with postprocess" $ do
    let cntx :: ROCContext Vec2 = clusterize [Vec2 0 0, Vec2 0.1 0] $ emptyROCContext cfg'
    assertBool "should be 2 clusters" $ length (nub $ rocPrototypes cntx) == 1
