from sortedcontainers import SortedDict, SortedList, SortedSet
from variant_candidate import VariantCandidate

class VariantCandidateCollection:
    def __init__(self):
        self.variants = SortedDict()
        
    def add_variant_candidate(self, variant_candidate):
        self.variants[variant_candidate.get_center()] = variant_candidate     
     
    def add_read_pair(self, read_pair, var_type):
        
        if len(self.variants) == 0:
            new_var_can = VariantCandidate(read_pair, var_type)
            self.add_variant_candidate(new_var_can)
        else:
            center = read_pair.get_center()
            ind = self.variants.bisect(center)
            if ind >= 1:
                var_key = self.variants.iloc[ind-1]
                var_can = self.variants[var_key]
    
                if var_can.var_type == var_type:
                    if var_can.belongs(read_pair):
                        var_can.add_read_pair(read_pair)
                        return
                else:
                    print "Multiple variants at same location"
                    
            if ind < len(self.variants):
                var_key = self.variants.iloc[ind]
                var_can = self.variants[var_key]
    
                if var_can.belongs(read_pair):
                    var_can.add_read_pair(read_pair)
                    return
 
            new_var_can = VariantCandidate(read_pair, var_type)
            self.add_variant_candidate(new_var_can)
    
    def merge(self):
        first = True
        for key in self.variants:
            if first:
                cur_var = self.variants[key]
                first = False
                continue
            candidate = self.variants[key]
            if cur_var.overlaps(candidate):
                cur_var.merge_with(candidate)
                del self.variants[key]
            else:
                cur_var = candidate
                
    def purge(self):
        keys_to_purge = []
        for key in iter(self.variants):
            if self.variants[key].get_inner_span() < 0 or self.variants[key].depth <= 1:
               keys_to_purge.append(key)
        for key in keys_to_purge:
            del self.variants[key]