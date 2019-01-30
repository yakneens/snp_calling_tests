from model.locus import Locus
import jsonpickle

l = Locus(pos=1211, ref=123454, alt=None, gl_ref=0, gl_het=0,
      gl_hom=0, history=[])

l2 = [l]

x = jsonpickle.encode(l)
y = jsonpickle.encode(l2)
print(str(jsonpickle.decode(x)))
print(jsonpickle.decode(y))