class VariantCandidate:
    def __init__(self, read_pair, var_type):
        self.supporting_reads = []
        self.supporting_reads.append(read_pair)
        
        self.left_pos = read_pair.get_left_pos()
        self.inner_left = read_pair.get_left_pos()
        self.right_pos = read_pair.get_right_pos()
        self.inner_right = read_pair.get_right_pos()
        self.var_type = var_type
        
        if self.left_pos > self.right_pos:
            print read_pair
            
        self.center = read_pair.get_center()
        
        self.width = read_pair.get_width()
        
        self.depth = 1
    
    def add_read_pair(self, read_pair):
        self.supporting_reads.append(read_pair)
        self.depth += 1
        self.left_pos = min(self.left_pos, read_pair.get_left_pos())
        self.inner_left = max(self.inner_left, read_pair.get_left_pos())
        
        self.right_pos = max(self.right_pos, read_pair.get_right_pos())
        self.inner_right = min(self.inner_right, read_pair.get_right_pos())
        self.center = (self.left_pos + self.right_pos) / 2
    
    def get_center(self):
        return self.center
    
    def get_width(self):
        return self.right_pos - self.left_pos
    
    def belongs(self, read_pair):
        the_center = read_pair.get_center()
        return the_center >= self.left_pos and the_center <= self.right_pos
        
    def overlaps(self, another):
        return (another.left_pos <= self.right_pos and another.left_pos >= self.left_pos) or (another.right_pos >= self.left_pos and another.right_pos <= self.right_pos)    
     
    def merge_with(self, another):
        for read_pair in another.supporting_reads:
            self.add_read_pair(read_pair)
            
    def inner_span_to_str(self):
        return "[inner-left:{} inner-right:{} span:{} depth:{}]".format(self.inner_left, self.inner_right, self.inner_right - self.inner_left, self.depth)   
            
    def get_inner_span(self):
        return self.inner_right - self.inner_left
    
    def __str__(self):
        ret = "{} ".format(self.var_type)
        for read in self.supporting_reads:
            ret += "[{},{}]\n".format(read.get_left_pos(), read.get_right_pos())
        return ret