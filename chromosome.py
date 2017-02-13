class chromosome:
    def __init__(self):
        self.windows = list()

class window:
    def __init__(self):
        self.variants = list()

class variant:
    def __init__(self, chrm, pos, ref, alt, nex, ntr, tex, ttr):
        self.chrm = chrm
        self.ref = ref
        self.alt = alt
        self.pos = pos
        self.nex = nex
        self.ntr = ntr
        self.tex = tex
        self.ttr = ttr
