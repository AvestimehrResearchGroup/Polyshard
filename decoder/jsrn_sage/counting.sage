"""Wrapper-algebra for counting no. of operations in fields"""
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.morphism import Morphism
from sage.structure.element import FieldElement
from sage.rings.ring import Field
# TODO:
# - Add explicit forwards to dir for tab completion
#
#
#
#

class ForwarderClass(object):
    # Forward all calls to this object for methods in self._allowed_forwards to self.obj

    def __init__(self):
        # By doing this instead of using __getattr__, we make sure that these
        # forwards override inherited versions of the methods
        for method in self._allowed_forwards:
            self.__dict__[method] = getattr(self.obj, method)

    def wrapped_obj(self):
        return obj

    def __hasattr__(self, name):
        #print "%s asked for attribute %s" % (self.__class__, name)
        return name in self._allowed_forwards

    def __dir__(self):
        return list(self._allowed_forwards)


class CountingFieldElement(ForwarderClass, FieldElement):
    _allowed_forwards = { '__str__', '__repr__', '_latex_', 'category', 'is_zero' }

    def __init__(self, CF, elem):
        self.CF = CF
        if isinstance(elem, int):
            elem = CF.obj(elem)
        self.obj = elem
        FieldElement.__init__(self, CF)
        ForwarderClass.__init__(self)

    def parent(self):
        return self.CF

    def _is_same_CF(self, other):
        return isinstance(other, self.__class__) and other.CF == self.CF

    def _is_over_CF(self, other):
        return hasattr(other, "base_ring") and other.base_ring() == self.CF
    
    def _add_(left, right):
        if not( left.obj.is_zero() or right.obj.is_zero() ):
            left.CF._additions += 1
        return left.CF.element_cache[left.obj + right.obj]

    def _iadd_(left, right):
        if not( left.obj.is_zero() or right.obj.is_zero() ):
            left.CF._additions += 1
        left.obj._iadd_(right.obj)
        
    def _sub_(left, right):
        if not( left.obj.is_zero() or right.obj.is_zero() ):
            left.CF._additions += 1
        return left.CF.element_cache[left.obj - right.obj]

    def _isub_(left, right):
        if not( left.obj.is_zero() or right.obj.is_zero() ):
            left.CF._additions += 1
        left.obj._isub_(right.obj)

    def _mul_(left, right):
        if not( left.obj.is_zero() or right.obj.is_zero() or left.obj.is_one() or right.obj.is_one() ):
            left.CF._multiplications += 1
        return left.CF.element_cache[left.obj * right.obj]
        
    def _imul_(left, right):
        if not( left.obj.is_zero() or right.obj.is_zero() or left.obj.is_one() or right.obj.is_one() ):
            left.CF._multiplications += 1
        left.obj._imul_(right.obj)

    def _div_(left, right):
        if not( left.obj.is_zero() or right.obj.is_zero() or left.obj.is_one() or right.obj.is_one() ):
            left.CF._multiplications += 1
        return left.CF.element_cache[left.obj / right.obj]
        
    def _idiv_(left, right):
        if not( left.obj.is_zero() or right.obj.is_zero() or left.obj.is_one() or right.obj.is_one() ):
            left.CF._multiplications += 1
        left.obj._idiv_(right.obj)

    def _latex_(self):
        return self.obj._latex_()

    def _repr_(self):
        return self.obj._repr_()

    def _str_(self):
        return self.obj._str_()

    def __nonzero__(self):
        return self.obj.__nonzero__()
        
    def __cmp__(left, right):
        return left.obj.__cmp__(right.obj)

    def __hash__(self):
        return self.obj.__hash__()

# Will look like a dictionary
class FakeCache():
    def __init__(self, CF):
        self.CF = CF
    def __getitem__(self, e):
        return CountingFieldElement(self.CF, e)
    def values(self):
        for e in self.CF.obj:
            yield CountingFieldElement(self.CF, e)

class CountingField(ForwarderClass, UniqueRepresentation, Field):
    _allowed_forwards = { 'characteristic', 'cardinality', 'is_finite' }

    _wrap_forwards = { 'next' }

    #TODO: roots

    Element = CountingFieldElement

    def __init__(self, F):
        Field.__init__(self, self)
        self.obj = F
        ForwarderClass.__init__(self)
        self._multiplications = 0
        self._additions = 0
        # parent stuff for coercion
        name = "CountingField"
        Parent.__init__(self)
        self._populate_coercion_lists_()
        self.rename(name)

        # Cache all elements
        if F.cardinality() < 10000:
            self.element_cache = { e : CountingFieldElement(self, e) for e in F }
        else:
            self.element_cache = FakeCache(self)

        # Other
        #self._zero_cache = self._element_constructor_(0)
        
    def base_ring(self):
        return self

    def additions(self):
        return self._additions

    def multiplications(self):
        return self._multiplications

    def reset_counters(self):
        ret = (self._additions, self._multiplications)
        self._additions = 0
        self._multiplications = 0
        return ret

    def pause_counters(self):
        self._counter_cache = (self._additions, self._multiplications)

    def resume_counters(self):
        (self._additions, self._multiplications) = self._counter_cache
        
    def __getattr__(self, name):
        # Apparently, these wrapper functions are not easily generated as the _allowed_forwards
        if name in self._wrap_forwards:
            def wrapper(*args, **nargs):
                return CountingFieldElement(self, getattr(self.obj, name)(*args, **nargs))
            return wrapper
        else:
            return super(CountingField, self).__getattr__(name)

    def __contains__(self, x):
        if isinstance(x, CountingFieldElement):
            return x.CF == self
        else:
            return False

    def __iter__(self):
        for e in self.element_cache.values():
            yield e

    def __getitem__(self,x):
        return PolynomialRing(self, x)

    # Coercion stuff
    def _coerce_map_from_(self, S):
        #print "Coercion counting field on: ", S
        return S==self.obj or S==CountingField or S==ZZ

    def _element_constructor_(self, e):
        e = self.obj._element_constructor_(e)
        return self.element_cache[e]
        
    def __eq__(self,other):
        return isinstance(other, CountingField) and self.obj == other.obj

    # Special stuff
    # Have to be overridden since Ring's default implementation are stupid
    def gen(self):
        return self.element_cache[self.obj.gen()]

    def random_element(self):
        return self.element_cache[self.obj.random_element()]

    def _repr_(self):
        return "Counting(%s)" % self.obj

    def __str__(self):
        return "Counting(%s)" % self.obj

    def _latex_(self):
        return "C^\#(%s)" % latex(self.obj)
