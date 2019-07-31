#Copyright Johan S. R. Nielsen
#This file is part of the codinglib for Sage.
#This library is free software released under the GPL v3 or any later version as
#published by the Free Software Foundation.

"""Abstract class for representing channels"""

from codinglib.util import *

class Channel(object):
    r"""
    Abstract class for representing channels.
    Channels "transmit" from one space to another, adding noise.

    FIELDS:

    - ``input_space`` - The space from which the input comes for this channel
    - ``output_space`` - The space in which this channel's output lands

    Subclasses should at least overwrite ``transmit_element''.
    """


    def __init__(self, input_space, output_space):
        self.input_space = input_space
        self.output_space = output_space

    def transmit_element(self, payload):
        """Transmit the payload being a single element in
        ``self.input_space''. Marginally faster than ``self.transmit''
        if one is transmitting a single element.
        Should raise TypeError if payload is not in ``self.input_space''.
        """
        raise NotImplementedError("Abstract base class Channel does not implement transmit_element.")


    def transmit(self, payload):
        """Transmit the given payload over the channel. The payload
        can be either one element from ``self.input_space'', or it can
        be a list of such"""
        if payload in self.input_space:
            return self.transmit_element(payload)
        else:
            # Assume it is an iterable of elements
            try:
                outs = []
                for p in payload:
                    outs.append(self.transmit_element(p))
                return outs
            except TypeError(e):
                raise TypeError("Channel can only transmit elements of input_space or iterables of such.")


class QarySymmetricChannel(Channel):
    r"""
    The q-ary symmetric channel over a finite alphabet

    FIELDS:

    - ``space`` - Some finite alphabet (with at least 2 elements), supporting a random_element method.
    - error_probability - Probability that a transmitted symbol will be modified by the channel
    """

    def __init__(self, space, error_probability):
        super(QarySymmetricChannel, self).__init__(space, space)
        self.space = space
        self.error_probability = error_probability

    def transmit_element(self, payload):
        if random() < self.error_probability:
            e = payload
            while e == payload:
                e = self.space.random_element() #TODO: should check for non-zero
            return e
        else:
            return payload


class BinaryAdditiveGaussianWhiteNoiseChannel(Channel):
    r"""
    The additive Gaussian white noise channel on a BSPK encoding (binary -1, 1).
    Input space is binary 0,1 and the conversion to BSPK is done by
    the channel before adding noise.
    The output space is the real numbers.
    """
    def __init__(self, sigma):
        super(BinaryAdditiveGaussianWhiteNoiseChannel, self).__init__(GF(2), RR)
        #TODO: Calculate sigma from SNR
        self._gaussian = RealDistribution('gaussian', sigma)

    def transmit_element(self, payload):
        bspk = 1 if payload else -1
        return bspk + self._gaussian.get_random_element()

    def decide_hard(self, received):
        """Return the most likely symbol given the received observation"""
        return self.input_space.zero() if received < 0 else self.input_space.one()

    def posteriori_probability(self, received, input):
        """Returns the a posteriori probability that ``input'' was sent given that ``received'' was observed"""
        bspk = 1 if input else -1
        dx  = received - bspk
        dnx = received + bspk
        gamma_x  = self._gaussian.distribution_function(dx)
        gamma_nx = self._gaussian.distribution_function(dnx)
        return gamma_x/(gamma_x+gamma_nx)

class BinaryExtSoftChannel(Channel):
    r"""
    A soft-information channel on a binary extension field, by
    representing a field element in binary and transmitting each
    symbol individually over a some soft channel.

    FIELDS

    - ``F``   - the binary extension field
    - ``bit_channel`` - the channel on which the bits are transmitted
    - ``m`` - number of bits in each symbol of F, i.e. $\log_2(|F|)$.
    """

    def __init__(self, F, bit_channel):
        assert(F.characteristic() == 2 and F.is_finite())
        self.m = F.degree()
        super(BinaryExtSoftChannel, self).__init__(F, VectorSpace(bit_channel.output_space, self.m))
        self.F = F
        self.bit_channel = bit_channel

    def transmit_element(self, payload):
        """Transmit the payload symbol by splitting it into bits and
        transmitting each bit over binary channel. Returns a list of
        observed bits."""
        return vector(self.bit_channel.transmit(payload._vector_()))

    def decide_hard(self, received):
        """Return the most likely symbol given the received observation"""
        hd_bits = [ self.bit_channel.decide_hard(ri) for ri in received ]
        return self.F(hd_bits)

    def posteriori_probability(self, received, input):
        """Returns the a posteriori probability that ``input'' was sent given that ``received'' was observed"""
        bits = input._vector_()
        return prod(self.bit_channel.posteriori_probability(r, b) for (r,b) in zip(received, bits))

    def posteriori_probability_table(self, received):
        """Calculate the entire table of all posteriori probabilities,
        given that ``received'' was observed. Returns a list of
        probabilities, corresponding to the symbols in the order of $F.list()$
        NOTE: This will take a long time and consume a lot of memory on very large fields."""
        zero_probs = list( self.bit_channel.posteriori_probability(r, 0) for r in received )
        one_probs =  list( 1-zp for zp in zero_probs )
        # TODO: Proceed in Gray code order
        #   gray_bin = (ii  >> 1) ^^ ii
        probs = []
        for a in self.F:
            bits = a._vector_()
            probs.append( prod( one_probs[i] if bits[i] else zero_probs[i] for i in range(self.m) ) )
        return probs
