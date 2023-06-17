from defines import SequenceAlignment

test = SequenceAlignment("jjjabcd", "abcdkk","global")

print(test)

print(test.reconstructGlobal())
