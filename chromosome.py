class chromosome:
    def __init__(self):
        self.windows = list()
        self.p_values = dict()
        self.emd_group0tex = dict()
        self.emd_group1tex = dict()
        self.emd_group0ttr = dict()
        self.emd_group1ttr = dict()
        self.bimodal={'tex': [], 'ttr': []}
        self.group=list()

    def join_groups(self, vaf):
        group = [[], []]
        for w in self.windows:
            if getattr(w, 'group_' + vaf) == 0:
                group[0] += w.variants
            elif getattr(w, 'group_' + vaf):
                group[1] += w.variants
        return group

    def join_windows(self, index1, index2):
        self.windows[index1].variants += self.windows[index2].variants
        self.windows[index2].variants = list()
        self.windows[index1].group.update([index1, index2])
        self.windows[index1].group.update(self.windows[index2].group)
        self.windows[index2].group = set()

    def getAllVariants(self):
        allVariants = list()
        for w in self.windows:
            allVariants += w.variants
        return allVariants

    def get_variants_number(self):
        for w in self.windows:
            w.var_num = len(w.variants)

class window:
    def __init__(self):
        self.variants = list()
        self.group = set()
        self.close_emd = None

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
