
from typing import Any


class SequenceAlignment:
    def __init__(self, sequenceOne,sequenceTwo, mismatchPenalty = -1, matchPenalty = 1, indelPenalty = -2):
        self.sequenceOne = sequenceOne
        self.sequenceTwo = sequenceTwo
        self.matchPenalty = matchPenalty
        self.mismatchPenalty = mismatchPenalty
        self.indelPenalty = indelPenalty
        self.pointers = None
        self.matrix = None
        self.reconstructedSequenceOne, self.reconstructedSequenceTwo = None,None

        
    def initMatrix(self):
        lenSequenceTwo, lenSequenceOne= len(self.sequenceTwo), len(self.sequenceOne)
        colLength, rowLength = lenSequenceOne+2, lenSequenceTwo+2
        matrix = [["-"]*rowLength for i in range((colLength))] 
        for i in range(lenSequenceOne):
            matrix[i+2][0] = self.sequenceOne[i]
            matrix[i+2][1] = (i+1)*self.indelPenalty
        for i in range(lenSequenceTwo):
            matrix[0][i+2] = self.sequenceTwo[i]
            matrix[1][i+2] = (i+1)*self.indelPenalty
        matrix[1][1] = 0
        self.matrix = matrix
    
    def __call__(self, *args: Any, **kwds: Any) -> Any:
        self.initMatrix()
        self.forwardPass()
        self.reconstruct()
        return self.matrix
    
    def __repr__(self) -> str:
        if not self.matrix:
            return '<%s.%s object at %s>' % (
                self.__class__.__module__,
                self.__class__.__name__,
                hex(id(self))
            )
        return ('\n'.join(' '.join(map(str, row)) for row in self.matrix) + '\n' + '\n'+
                 self.reconstructedSequenceOne + "\n" + self.reconstructedSequenceTwo)

    


class GlobalAlignment(SequenceAlignment):
    def __init__(self, sequenceOne,sequenceTwo, mismatchPenalty = -1, matchPenalty = 1, indelPenalty = -2):
        super().__init__(sequenceOne,sequenceTwo, mismatchPenalty, matchPenalty, indelPenalty)


    def forwardPass(self):
        lenSequenceTwo, lenSequenceOne= len(self.sequenceTwo), len(self.sequenceOne)
        colLength, rowLength = lenSequenceOne+2, lenSequenceTwo+2
        pointers = dict()

        for i in range(1,colLength):
            for j in range(1,rowLength):
                if self.matrix[i][j] != "-":
                    continue
                
                indelPenalty = max(self.matrix[i-1][j], self.matrix[i][j-1]) + self.indelPenalty 
                misMatchPenalty = self.matrix[i-1][j-1] + self.matchPenalty if self.matrix[i][0] == self.matrix[0][j] else self.matrix[i-1][j-1] + self.mismatchPenalty
                self.matrix[i][j] = max(indelPenalty, misMatchPenalty)

                if self.matrix[i][j] == self.matrix[i-1][j]+ self.indelPenalty :
                    pointers[(i,j)] = [i-1,j]
                elif self.matrix[i][j] == self.matrix[i][j-1]+ self.indelPenalty :
                    pointers[(i,j)] = [i,j-1]
                else:
                    pointers[(i,j)] = [i-1,j-1]
        self.pointers = pointers


    def reconstruct(self):
        reconstructedOne = ""
        reconstructedTwo = ""
        rowIdx = len(self.matrix)-1
        colIdx = len(self.matrix[0])-1
        while True:
            if (rowIdx == 1) and (colIdx == 1):
                break
            rowIdxAbove = rowIdx-1
            colIdxLeft = colIdx-1

            if rowIdxAbove == 0:
                reconstructedTwo = self.matrix[0][colIdx] + reconstructedTwo
                reconstructedOne = "-" + reconstructedOne
                colIdx -=1
                continue
            if colIdxLeft == 0:
                reconstructedOne = self.matrix[rowIdx][0] + reconstructedOne
                reconstructedTwo = "-" + reconstructedTwo
                rowIdx -=1
                continue

            newRowIdx, newColIdk = self.pointers[(rowIdx,colIdx)]
            
            if newRowIdx != rowIdx and newColIdk == colIdx: 
                reconstructedOne = self.matrix[rowIdx][0] + reconstructedOne
                reconstructedTwo = "-" + reconstructedTwo
            elif newColIdk != colIdx and newRowIdx == rowIdx:
                reconstructedOne = "-" + reconstructedOne
                reconstructedTwo = self.matrix[0][colIdx] + reconstructedTwo
            else:
                reconstructedOne = self.matrix[rowIdx][0] + reconstructedOne
                reconstructedTwo = self.matrix[0][colIdx] + reconstructedTwo

            rowIdx = newRowIdx
            colIdx = newColIdk

        self.reconstructedSequenceOne, self.reconstructedSequenceTwo = reconstructedOne, reconstructedTwo
        
class LocalAlignment(SequenceAlignment):
    def __init__(self, sequenceOne,sequenceTwo, mismatchPenalty = -1, matchPenalty = 1, indelPenalty = -2):
        super().__init__(sequenceOne,sequenceTwo, mismatchPenalty, matchPenalty, indelPenalty)
        
    
    def forwardPass(self):
        ### height corresponds to sequenceOne, width is sequenceTwo
        ### [ [s2],// [s2],// [s2] ]
        lenSequenceTwo, lenSequenceOne= len(self.sequenceTwo), len(self.sequenceOne)
        colLength, rowLength = lenSequenceOne+2, lenSequenceTwo+2
        pointers = dict()
        for i in range(1,colLength):
            for j in range(1,rowLength):
                if self.matrix[i][j] != "-":
                    continue
                
                indelPenalty = max(self.matrix[i-1][j], self.matrix[i][j-1]) + self.indelPenalty 
                misMatchPenalty = self.matrix[i-1][j-1] + self.matchPenalty if self.matrix[i][0] == self.matrix[0][j] else self.matrix[i-1][j-1] + self.mismatchPenalty
                self.matrix[i][j] = max(indelPenalty, misMatchPenalty, 0)

                if self.matrix[i][j] == 0:
                    continue
                elif self.matrix[i][j] == self.matrix[i-1][j]+ self.indelPenalty :
                    pointers[(i,j)] = [i-1,j]
                elif self.matrix[i][j] == self.matrix[i][j-1]+ self.indelPenalty :
                    pointers[(i,j)] = [i,j-1]
                else:
                    pointers[(i,j)] = [i-1,j-1]
        self.pointers = pointers



    def reconstruct(self):
        reconstructedOne=""
        reconstructedTwo=""
        rowLength = len(self.matrix[0])
        colLength = len(self.matrix)
        maxEntry = 0
        iMax, jMax = 0 , 0
        for i in range(1,colLength):
            for j in range(1,rowLength):
                if self.matrix[i][j] > maxEntry:
                    maxEntry = self.matrix[i][j]
                    iMax, jMax = i, j

        rowIdx,colIdx = iMax, jMax

        while True:
            if self.matrix[rowIdx][colIdx] == 0:
                break
            rowIdxAbove = rowIdx-1
            colIdxLeft = colIdx-1

            if rowIdxAbove == 0:
                reconstructedTwo = self.matrix[0][colIdx] + reconstructedTwo
                reconstructedOne = "-" + reconstructedOne
                colIdx -=1
                continue
            if colIdxLeft == 0:
                reconstructedOne = self.matrix[rowIdx][0] + reconstructedOne
                reconstructedTwo = "-" + reconstructedTwo
                rowIdx -=1
                continue

            newRowIdx, newColIdk = self.pointers[(rowIdx,colIdx)]
            
            if newRowIdx != rowIdx and newColIdk == colIdx: 
                reconstructedOne = self.matrix[rowIdx][0] + reconstructedOne
                reconstructedTwo = "-" + reconstructedTwo
            elif newColIdk != colIdx and newRowIdx == rowIdx:
                reconstructedOne = "-" + reconstructedOne
                reconstructedTwo = self.matrix[0][colIdx] + reconstructedTwo
            else:
                reconstructedOne = self.matrix[rowIdx][0] + reconstructedOne
                reconstructedTwo = self.matrix[0][colIdx] + reconstructedTwo

            rowIdx = newRowIdx
            colIdx = newColIdk

        self.reconstructedSequenceOne, self.reconstructedSequenceTwo = reconstructedOne, reconstructedTwo
        
