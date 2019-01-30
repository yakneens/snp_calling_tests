class ReadPair:
    def __init__(self, read_id, first, second):
        
        self.first = first
        self.second = second
            
        self.read_id = read_id
     
    def get_center(self):
        return (self.first.reference_start + self.second.reference_start + self.second.query_alignment_length) / 2
    
    def get_width(self):
        return self.second.reference_start + self.second.query_alignment_length - self.first.reference_start
           
    def get_left_pos(self):
        return self.first.reference_start
    
    def get_right_pos(self):       
        return self.second.reference_start + self.second.query_alignment_length
    
    def get_insert_size(self):
        return self.second.reference_start - self.first.reference_start
    
    def __str__(self):
        return "{} {} {} {}".format(self.read_id, self.get_center(), self.first, self.second)
    
    def was_reversed(self):
        return self.was_reversed
    
    def reads_overlap(self):
        return True if self.first.mpos < self.first.aend else False