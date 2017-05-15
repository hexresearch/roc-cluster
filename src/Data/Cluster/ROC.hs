module Data.Cluster.ROC(
  -- * Algorithm configuration
    ROCConfig
  , rocThreshold
  , rocMaxClusters
  , defaultROCConfig
  -- * Cluster definition
  , Prototype
  , prototypeValue
  , prototypeWeight
  -- * API
  , ClusterSpace(..)
  , ROCContext
  , emptyROCContext
  , loadROCContext
  , rocPrototypes
  , clusterize
  ) where

import Data.Data
import Data.Monoid
import Data.Ord
import Data.Vector (Vector)
import GHC.Generics

import qualified Data.Foldable as F
import qualified Data.Vector as V

-- | Configuration of ROC clusterization
data ROCConfig = ROCConfig {
  -- | If weight of prototype is less than the value, it is removed at final
  -- step.
  rocThreshold   :: !Double
  -- | Maximum count of clusters, could be less
, rocMaxClusters :: !Int
} deriving (Generic, Data)

-- | Default configuration:
-- @
-- ROCConfig {
--   rocThreshold = 0
-- , rocMaxClusters = 10
-- }
-- @
defaultROCConfig :: ROCConfig
defaultROCConfig = ROCConfig {
    rocThreshold = 0
  , rocMaxClusters = 10
  }

-- | Operations that value has to support to use in ROC clusterisation
class ClusterSpace a where
  -- | Zero point in space
  pointZero   :: a
  -- | Addition of vectors in space
  pointAdd    :: a -> a -> a
  -- | Scaling by a scalar
  pointScale  :: Double -> a -> a
  -- | Kernel function
  pointKernel :: a -> a -> Double
  -- | Square of distance between of points (defined via kernel) and exposed
  -- only for possible optimizations as for Gaussian kernel (2 - 2 * pointKernel x y)
  pointDistanceSquared :: a -> a -> Double
  pointDistanceSquared x y = pointKernel x x - 2 * pointKernel x y + pointKernel y y
  {-# INLINE pointDistanceSquared #-}

-- | Cluster information
data Prototype a = Prototype {
  prototypeValue  :: !a
, prototypeWeight :: !Double
} deriving (Generic, Functor)

instance ClusterSpace a => Monoid (Prototype a) where
  mempty = Prototype pointZero 0
  mappend p1 p2 = Prototype pos w
    where
      w = prototypeWeight p1 + prototypeWeight p2
      pos = (1/w) `pointScale` ((prototypeWeight p1 `pointScale` prototypeValue p1) `pointAdd` (prototypeWeight p2 `pointScale` prototypeValue p2))
  {-# INLINE mempty #-}
  {-# INLINE mappend #-}

-- | Internal context of algorithm
data ROCContext a = ROCContext {
  cntxPrototypes :: !(Vector (Prototype a))
, cntxConfig     :: !ROCConfig
} deriving (Generic, Functor)

-- | Create new context for clusterization from scratch
emptyROCContext :: ROCConfig -> ROCContext a
emptyROCContext cfg = ROCContext {
    cntxPrototypes = mempty
  , cntxConfig     = cfg
  }

-- | Load context from set of prototypes
loadROCContext :: Foldable f => ROCConfig -> f (Prototype a) -> ROCContext a
loadROCContext cfg ps = (emptyROCContext cfg) { cntxPrototypes = V.fromList . F.toList $ ps }

-- | Get collection of prototypes from ROC context
rocPrototypes :: ROCContext a -> [Prototype a]
rocPrototypes = F.toList . cntxPrototypes

-- | Perform clusterization of next part of data
clusterize :: forall a f . (ClusterSpace a, Foldable f)
  => f a -- ^ Set of data that need to be added to clusters
  -> ROCContext a -- ^ Context with current prototypes
  -> ROCContext a -- ^ Updated context
clusterize xs cntx0 = clusterizePostprocess addAll
  where
    addAll = F.foldl' (flip clusterizeAddMerge) cntx0 xs

-- | Cluster a single value (step 2-6 in original paper). Moves existing clusters,
-- creates new clusters and merges close clusters.
clusterizeAddMerge :: forall a . (ClusterSpace a)
  => a -- ^ Single point
  -> ROCContext a -- ^ Context with current prototypes
  -> ROCContext a -- ^ Updated context
clusterizeAddMerge x cntx = clusterizeNewPrototype x $ if n >= nmax then clusterizeMerge cntx' else cntx'
  where
    cntx' = clusterizeSingle x cntx
    n = V.length . cntxPrototypes $ cntx'
    nmax = rocMaxClusters . cntxConfig $ cntx'
{-# INLINE clusterizeAddMerge #-}

-- | Cluster a single value (step 2 in original paper). This step updates only existing
-- clusters.
clusterizeSingle :: forall a . (ClusterSpace a)
  => a -- ^ Single point
  -> ROCContext a -- ^ Context with current prototypes
  -> ROCContext a -- ^ Updated context
clusterizeSingle x ctx@ROCContext{..}
  | V.null cntxPrototypes = ctx
  | otherwise             = ctx { cntxPrototypes = cntxPrototypes V.// [(winnerIndex, winner')] }
    where
    winnerIndex = V.minIndex . fmap (pointDistanceSquared x . prototypeValue) $ cntxPrototypes
    winner = cntxPrototypes V.! winnerIndex
    winner' = let Prototype{..} = winner in Prototype
      (prototypeValue `pointAdd` ( (1 / prototypeWeight) `pointScale` (x `pointAdd` pointScale (-1) prototypeValue) ))
      (prototypeWeight + pointKernel x prototypeValue)
{-# INLINE clusterizeSingle #-}

-- | Merge the most closest clusters (step 4 in original paper).
clusterizeMerge :: forall a . (ClusterSpace a)
  => ROCContext a -- ^ Context with current prototypes
  -> ROCContext a -- ^ Updated context
clusterizeMerge ctx@ROCContext{..}
  | V.length cntxPrototypes <= 1 = ctx
  | otherwise = ctx { cntxPrototypes = cntxPrototypes' }
  where
    -- find two prototypes that have minimum distance (warning, Vector monad!)
    (minxi, minyi, _) = V.minimumBy (comparing $ \(_, _, a) -> a) $ do
      (xi, xv) <- V.indexed cntxPrototypes
      (yi, yv) <- V.take xi $ V.indexed cntxPrototypes
      pure (xi, yi, prototypeValue yv `pointDistanceSquared` prototypeValue xv)
    x = cntxPrototypes V.! minxi
    y = cntxPrototypes V.! minyi
    x' = x <> y
    removeAt i v = V.slice 0 i v <> V.slice (i+1) (V.length v - i - 1) v
    cntxPrototypes' = removeAt minyi $ cntxPrototypes V.// [(minxi, x')]
{-# INLINE clusterizeMerge #-}

-- | Form a new prototype from single point (step 5 in original paper)
clusterizeNewPrototype :: forall a . (ClusterSpace a)
  => a -- ^ Point
  -> ROCContext a -- ^ Context with current prototypes
  -> ROCContext a -- ^ Updated context
clusterizeNewPrototype a ctx@ROCContext{..} = ctx { cntxPrototypes = cntxPrototypes `V.snoc` newProto }
  where
    newProto = Prototype a 0
{-# INLINE clusterizeNewPrototype #-}

-- | Remove clusters that have negligible weights (step 6 in original paper)
clusterizePostprocess :: forall a . (ClusterSpace a)
  => ROCContext a -- ^ Context with current prototypes
  -> ROCContext a -- ^ Updated context
clusterizePostprocess ctx@ROCContext{..} = ctx { cntxPrototypes = V.filter isValuable cntxPrototypes }
  where
    threshold = rocThreshold cntxConfig
    isValuable p = prototypeWeight p > threshold
{-# INLINE clusterizePostprocess #-}
