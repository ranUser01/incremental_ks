from random import random

class Treap:
  def __init__(self, key, value = 0):
    self.key = key
    self.value = value
    self.priority = random()
    self.size = 1
    self.height = 1
    self.lazy = 0
    self.max_value = value
    self.min_value = value
    self.left = None
    self.right = None

  @staticmethod
  def SumAll(node, value):
    if node is None:
      return
    node.value += value
    node.max_value += value
    node.min_value += value
    node.lazy += value
  
  @classmethod
  def Unlazy(cls, node):
    cls.SumAll(node.left, node.lazy)
    cls.SumAll(node.right, node.lazy)
    node.lazy = 0

  @classmethod
  def Update(cls, node):
    if node is None:
      return
    cls.Unlazy(node)
    node.size = 1
    node.height = 0
    node.max_value = node.value
    node.min_value = node.value

    if node.left is not None:
      node.size += node.left.size
      node.height = node.left.height
      node.max_value = max(node.max_value, node.left.max_value)
      node.min_value = min(node.min_value, node.left.min_value)
    
    if node.right is not None:
      node.size += node.right.size
      node.height = max(node.height, node.right.height)
      node.max_value = max(node.max_value, node.right.max_value)
      node.min_value = min(node.min_value, node.right.min_value)
    
    node.height += 1
  
  @classmethod
  def SplitKeepRight(cls, node, key):
    if node is None:
      return None, None

    left, right = None, None
    
    cls.Unlazy(node)

    if key <= node.key:
      left, node.left = cls.SplitKeepRight(node.left, key)
      right = node
    else:
      node.right, right = cls.SplitKeepRight(node.right, key)
      left = node
    
    cls.Update(left)
    cls.Update(right)

    return left, right

  @classmethod
  def Merge(cls, left, right):
    if left is None:
      return right
    if right is None:
      return left
    
    node = None

    if left.priority > right.priority:
      cls.Unlazy(left)
      left.right = cls.Merge(left.right, right)
      node = left
    else:
      cls.Unlazy(right)
      right.left = cls.Merge(left, right.left)
      node = right
    
    cls.Update(node)
    return node

  @classmethod
  def SplitSmallest(cls, node):
    if node is None:
      return None, None

    left, right = None, None
    
    cls.Unlazy(node)

    if node.left is not None:
      left, node.left = cls.SplitSmallest(node.left)
      right = node
    else:
      right = node.right
      node.right = None
      left = node
    
    cls.Update(left)
    cls.Update(right)

    return left, right

  @classmethod
  def SplitGreatest(cls, node):
    if node is None:
      return None, None
    
    cls.Unlazy(node)

    if node.right is not None:
      node.right, right = cls.SplitGreatest(node.right)
      left = node
    else:
      left = node.left
      node.left = None
      right = node
    
    cls.Update(left)
    cls.Update(right)

    return left, right

  @staticmethod
  def Size(node):
    return 0 if node is None else node.size
  
  @staticmethod
  def Height(node):
    return 0 if node is None else node.height

  @classmethod
  def _ToList(cls, node, extractor, _list = None):
    if _list is None:
      _list = []
    if node is None:
      return _list
    cls.Unlazy(node)
    cls._ToList(node.left, extractor, _list)
    _list.append(extractor(node))
    cls._ToList(node.right, extractor, _list)
    return _list
  
  @classmethod
  def KeysToList(cls, node, _list = None):
    extractor = lambda x: x.key
    return cls._ToList(node, extractor, _list)

  @classmethod
  def ValuesToList(cls, node, _list = None):
    extractor = lambda x: x.value
    return cls._ToList(node, extractor, _list)

from math import log

class IKS:
  def __init__(self):
    self.treap = None
    self.n = [0, 0]

  @staticmethod
  def KSThresholdForPValue(pvalue, N):
    '''Threshold for KS Test given a p-value
    Args:
      pval (float): p-value.
      N (int): the size of the samples.

    Returns:
      Threshold t to compare groups 0 and 1. The null-hypothesis is discarded if KS() > t.
    '''
    ca = (-0.5 * log(pvalue)) ** 0.5
    return ca * (2.0 * N / N ** 2)

  @staticmethod
  def CAForPValue(pvalue):
    '''ca for KS Test given a p-value
    Args:
      pval (float): p-value.

    Returns:
      Threshold the "ca" that can be used to compute a threshold for KS().
    '''
    return (-0.5 * log(pvalue)) ** 0.5
  
  def KS(self):
    '''Kolmogorov-Smirnov statistic. Both groups must have the same number of observations.

    Returns:
      The KS statistic D.
    '''
    assert(self.n[0] == self.n[1])
    N = self.n[0]
    if N == 0:
      return 0
    return max(self.treap.max_value, -self.treap.min_value) / N
  
  def Kuiper(self):
    '''Kuiper statistic. Both groups must have the same number of observations.

    Returns:
      The Kuiper statistic.
    '''
    assert(self.n[0] == self.n[1])
    N = self.n[0]
    if N == 0:
      return 0
    return (self.treap.max_value - self.treap.min_value) / N
  
  def Add(self, obs, group):
    '''Insert new observation into one of the groups.

    Args:
      obs: the value of the obseration. Tip: a tuple (actual value, random value) is recommended when there is overlap between groups or if values are not guaranteed to be mostly unique.
      group (int): which group the observation belongs to. Must be either 0 or 1.
    '''
    group = 0 if group == 2 else group
    assert(group == 0 or group == 1)
    key = (obs, group)

    self.n[group] += 1
    left, left_g, right, val = None, None, None, None

    left, right = Treap.SplitKeepRight(self.treap, key)

    left, left_g = Treap.SplitGreatest(left)
    val = 0 if left_g is None else left_g.value
    left = Treap.Merge(left, left_g)

    right = Treap.Merge(Treap(key, val), right)

    Treap.SumAll(right, 1 if group == 0 else -1)

    self.treap = Treap.Merge(left, right)

  def Remove(self, obs, group):
    '''Remove observation from one of the groups.

    Args:
      obs: the value of the obseration. Must be identical to a previously inserted observation (including the random element of a tuple, if this was the case).
      group (int): which group the observation belongs to. Must be either 0 or 1.
    '''
    group = 0 if group == 2 else group
    assert(group == 0 or group == 1)
    key = (obs, group)

    self.n[group] -= 1

    left, right, right_l = None, None, None

    left, right = Treap.SplitKeepRight(self.treap, key)
    right_l, right = Treap.SplitSmallest(right)

    if right_l is not None and right_l.key == key:
      Treap.SumAll(right, -1 if group == 0 else 1)
    else:
      right = Treap.Merge(right_l, right)
    
    self.treap = Treap.Merge(left, right)
  
  def Test(self, ca = 1.95):
    '''Test whether the reference and sliding window follow the different probability distributions according to KS Test.

    Args:
      ca: ca is a parameter used to calculate the threshold for the Kolmogorov-Smirnov statistic. The default value corresponds to a p-value of 0.001. Use IKS.CAForPValue to obtain an appropriate ca.

    Returns:
      True if we **reject** the null-hypothesis that states that both windows have the same distribution. In other words, we can consider that the windows have now different distributions.
    '''
    ca = ca or 1.95
    n = self.n[0]
    return self.KS() > ca * (2 * n / n ** 2) ** 0.5

IKS.AddObservation = IKS.Add
IKS.RemoveObservation = IKS.Remove

from collections import deque

class IKSSW:
  def __init__(self, values):
    '''Incremental Kolmogorov-Smirnov Sliding Window. This class assumes that one window is fixed (reference window) and another slides over a stream of data. The reference window can be updated to be the same as the current sliding window.

    Args:
      values: initial values for the reference and sliding windows.
    '''
    self.iks = IKS()
    self.sw = deque()
    self.reference = [(x, random()) for x in values]

    for val in self.reference:
      self.iks.AddObservation(val, 1)

    for val in values:
      wrnd = (val, random())
      self.sw.append(wrnd)
      self.iks.AddObservation(wrnd, 2)

  def Increment(self, value):
    '''Remove the oldest observation from the sliding window and replace it with a given value.
    
    Args:
      value: the new observation.
    '''
    self.iks.RemoveObservation(self.sw.popleft(), 2)
    wrnd = (value, random())
    self.iks.AddObservation(wrnd, 2)
    self.sw.append(wrnd)

  __call__ = Increment

  def Kuiper(self):
    '''Kuiper statistic. Both groups must have the same number of observations.

    Returns:
      The Kuiper statistic.
    '''
    return self.iks.Kuiper()

  def KS(self):
    '''Kolmogorov-Smirnov statistic. Both groups must have the same number of observations.

    Returns:
      The KS statistic D.
    '''
    return self.iks.KS()

  def Update(self):
    '''Updates the IKSSW. The reference window becomes the sliding window.
    '''
    for val in self.reference:
      self.iks.Remove(val, 1)

    self.reference.clear()
    for x in self.sw:
      self.reference.append((x[0], random()))

    for val in self.reference:
      self.iks.Add(val, 1)
  
  def Test(self, ca = 1.95):
    '''Test whether the reference and sliding window follow the different probability distributions according to KS Test.

    Args:
      ca: ca is a parameter used to calculate the threshold for the Kolmogorov-Smirnov statistic. The default value corresponds to a p-value of 0.001. Use IKS.CAForPValue to obtain an appropriate ca.

    Returns:
      True if we **reject** the null-hypothesis that states that both windows have the same distribution. In other words, we can consider that the windows have now different distributions.
    '''
    return self.iks.Test(ca)

if __name__ == "__main__":
  v = [random() for x in range(10)]
  ikssw = IKSSW(v)
  print(ikssw.KS(), ikssw.Kuiper(), ikssw.Test())
  for i in range(10):
    ikssw(random())
  print(ikssw.KS(), ikssw.Kuiper(), ikssw.Test())
  ikssw.Update()
  print(ikssw.KS(), ikssw.Kuiper(), ikssw.Test())
