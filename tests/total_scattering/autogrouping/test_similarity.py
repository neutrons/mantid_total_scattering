# test it
import unittest
import numpy
import total_scattering.autogrouping.similarity as similarity

class test_similarity( unittest.TestCase ):

  # return r'th component of a triangle weight function extending L indices in each direction
  # from the "origin". Add 1 to L to account for the zero-offset ( r == 0 ) weight value
  def weight_triangular( self, r, L ):

    return ( L - abs( r ) + 1 ) / ( L + 1 )
  
  '''
  Weighted cross correlation of vectors within a small window:
  sum_over_r( weight[ r ] * sum_over_x ( f( x ) * g( x + r ) ) ) for r in
  union([ 0, L ] and [ 0, -L ] ), inclusive. By default, the weight function
  is the triangular weight function defined above. Tested with a known analytic
  solution calculated using numpy functions
  '''
  def test_weighted_cross_corr( self ):

    # Arbitrary vectors
    f = [ 1, 3, 5, 7, 11, 13, 15, 17 ]
    g = [ 2, 4, 6, 8, 10, 12, 14, 16 ]

    # Ensure that each component of f and g contributes to the result 
    L = len( f ) - 1

    weights = [ self.weight_triangular( r, L ) for r in range( -L, L + 1 ) ]

    # Numpy expects the vectors in the reverse order
    ncc = numpy.correlate( g, f, mode = 'full' )

    # Element-wise multiply to apply the weights:
    weighted_ncc = weights * ncc

    # Sum to get non-normalized similarity measure:
    gold = numpy.sum( weighted_ncc )

    # Calculate using the method under test:
    metric = similarity.similarity_metric( L )

    test = metric.triangular_weighted_cross_correlation( f, g )

    self.assertEqual( gold, test )

  '''
  Test the normalized similarity meausure built from the weighted cross correlation.
  It should be:
  commutative: similarity( f, g ) = similarity( g, f )
  unity for auto-similarity: similarity( f, f ) = 1
  indicate similarity: similarity( f, h ) > similarity( g, h ) ---> f and h are more similar
  '''
  def test_de_gelder_similarity( self ):

    # instantiate similarity object
    metric = similarity.similarity_metric()

    # the metric should return unity for perfect correlation
    f = [ 1, 2, 3, 4, 5, 6, 7, 8 ]

    metric_value = metric.de_gelder_similarity( f, f )

    self.assertEqual( metric_value, 1 )

    # since it is a normalized metric, it should be invariant to scaling.
    # it should also be commutative
    f = [ 1, 3, 5, 7, 11, 13 ]

    g = [ 2, 4, 6, 8, 10, 12 ]

    scale_factor = 5

    self.assertEqual( metric.de_gelder_similarity( scale_factor * f, scale_factor * g ), \
                      metric.de_gelder_similarity( scale_factor * g, scale_factor * f ) )

    # The metric should indicate similarity. This vector should be more similar to f than g
    h = [ x - 0.1 for x in f ]
    fh = metric.de_gelder_similarity( f, h )
    fg = metric.de_gelder_similarity( f, g )

    self.assertGreater( fh, fg )

  '''
  test equality via a known analytic solution
  '''
  def test_pointwise_squared_difference_similarity( self ):

    # instantiate similarity object
    metric = similarity.similarity_metric()

    diff = 2

    f = [ 1, 2, 3, 4, 5, 6, 7, 8 ]
    g = [ x + diff for x in f ]

    squared_diff_sum = metric.pointwise_squared_difference_similarity( f, g )

    gold = 1 - len( f ) * diff ** 2

    self.assertEqual( squared_diff_sum, gold )

if __name__ == '__main__':
  unittest.main()
