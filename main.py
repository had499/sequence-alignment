from sequence import GlobalAlignment, LocalAlignment

## Define the alignment object using input sequences
local_alignment = LocalAlignment("ATCGGCTAGCTAGGCCAAATCGAC","AGGTCGACAGGTCGACAGGTCGAC")

global_alignment = GlobalAlignment("ATCGGCTAGCTAGGCCAAATCGAC","AGGTCGACAGGTCGACAGGTCGAC")

local_alignment()
global_alignment()


print(local_alignment)
print(global_alignment)

