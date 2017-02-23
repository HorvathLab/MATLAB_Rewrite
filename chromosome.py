class chromosome:
    def __init__(self):
        self.windows = list()
        self.p_values = dict()
        self.emd_group0tex = dict()
        self.emd_group1tex = dict()
        self.emd_group0ttr = dict()
        self.emd_group1ttr = dict()

    def join_groups(self, vaf):
        group = [[], []]
        for w in self.windows:
            if getattr(w, 'group_' + vaf) == 0:
                group[0] += w.variants
            elif getattr(w, 'group_' + vaf):
                group[1] += w.variants
        return group

    def getAllVariants(self):
        allVariants = list()
        for w in self.windows():
            allVariants += w.variants
        return allVariants

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
