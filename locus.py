import sys
from math import log10


class Locus:
    MIN_GL = sys.float_info.min

    def __init__(self, pos, ref, alt, gl_ref, gl_het, gl_hom, history):
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.gl_ref = gl_ref
        self.gl_het = gl_het
        self.gl_hom = gl_hom
        self.history = history
        # Read Depth
        self.dp = 0
        # Reference observation count
        self.ro = 0
        # Sum of reference observation qualities
        self.qr = 0
        self.ao = 0
        # Sum of alternate observation qualities
        self.qa = 0

    def add_to_history(self, new_val):
        self.history.append(new_val)

    def get_gls(self):
        return (self.gl_ref, self.gl_het, self.gl_hom)

    def get_log_gls(self):
        return (log10(self.gl_ref), log10(self.gl_het), log10(self.gl_hom))

    def get_pl(self):
        vals = (int(to_phred(self.gl_ref)), int(to_phred(self.gl_het)), int(to_phred(self.gl_hom)))
        return map(lambda x: x - min(vals), vals)

    def get_gq(self):
        pl = self.get_pl()
        if (pl[1] == 0):
            return pl[2]
        elif (pl[2] == 0):
            return pl[1]
        else:
            print "Something went wrong in PL calculation"
            exit(1)

    def get_norm_gls(self):
        norm_const = self.gl_ref + self.gl_het + self.gl_hom
        return (self.gl_ref / norm_const, self.gl_het / norm_const, self.gl_hom / norm_const)

    def update_gls_to_min_if_zero(self):
        self.gl_hom = min_gl_if_zero(self.gl_hom)
        self.gl_het = min_gl_if_zero(self.gl_het)
        self.gl_ref = min_gl_if_zero(self.gl_ref)

def to_phred(x):
    return -10*(log10(x))

def from_phred(x):
    return pow(10, -x/10)

def min_gl_if_zero(val):
    if val == 0:
        return Locus.MIN_GL
    else:
        return val