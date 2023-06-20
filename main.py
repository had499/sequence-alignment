from defines import SequenceAlignment

test = SequenceAlignment("what the hell is going on", "i dont what the hell is going on")

print(test.reconstructGlobal())
print(test)
