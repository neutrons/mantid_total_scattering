# code it
# Bundle all of this into a class, rename it something better
import numpy

'''
This file calculates similarity metrics for numpy vectors.
It is build from cross correlations
'''


class similarity_metric:

    def __init__(self, L=5):

        self.L = L
        self.weights = [self.weight_triangular(r) for r in range(0, self.L + 1)]

    '''
    symmetric weight function
    '''

    def weight_triangular(self, r):

        return (self.L - abs(r) + 1) / (self.L + 1)

    '''
    f, g ---> vectors to cross correlate.

    L ---> From their initial aligned positions, vectors will be cross
           correlated at offsets in union([ 0, L ] and [ 0, -L ] ), inclusive.

    weights ---> weight[ i ] to be applied to cross_correlation[ i ]
    for i in [ -L, L ] inclusive

    returns ---> weighted cross correlation of the vectors:
                 sum_over_r( weight[ r ] * sum_over_x ( f( x ) * g( x + r ) ) )
    '''

    def triangular_weighted_cross_correlation(self, f, g):

        if (len(f) != len(g)):
            raise RuntimeError(
                "weighted_cross_correlation expects equal length vectors")

        # Calculate sum at r = 0
        weighted_sum = numpy.dot(f, g)

        for r in range(1, self.L + 1):
            positive_sub_sum = numpy.dot(f[: -r], g[r:])

            negative_sub_sum = numpy.dot(f[r:], g[: -r])

            weighted_sum += (
                    self.weights[r] * (positive_sub_sum + negative_sub_sum))

        return weighted_sum

    '''
    normalize the triangular weighted cross correlation to provide a
    pairwise similarity metric

    abs( metric ) <= 1
    '''

    def de_gelder_similarity(self, f, g):

        weighted_cross_correlation = \
            self.triangular_weighted_cross_correlation(f, g)

        g_norm_component = \
            self.triangular_weighted_cross_correlation(g, g)

        f_norm_component = \
            self.triangular_weighted_cross_correlation(f, f)

        normalization = numpy.sqrt(g_norm_component * f_norm_component)

        return weighted_cross_correlation / normalization

    '''
    a very simple similarity metric
    '''

    def pointwise_squared_difference_similarity(self, f, g):

        if (len(f) != len(g)):
            raise RuntimeError(
                "weighted_cross_correlation expects equal length vectors")

        return 1 - numpy.sum(numpy.square(numpy.subtract(f, g)))
