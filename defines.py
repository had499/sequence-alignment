

class SequenceAlignment:
    def __init__(self,sequenceOne,sequenceTwo, mismatchPenalty = -1, matchPenalty = 1, indelPenalty = -2):
        self.matchPenalty = matchPenalty
        self.mismatchPenalty = mismatchPenalty
        self.indelPenalty = indelPenalty
        self.globalPointers = None
        self.localPointers = None
        self.matrixGlobal = self.initMatrixGlobal(sequenceOne,sequenceTwo)
        self.matrixLocal = self.initMatrixLocal(sequenceOne,sequenceTwo)
        self.reconstructedSequenceOneGlobal, self.reconstructedSequenceTwoGlobal = self.reconstructGlobal() 
        self.reconstructedSequenceOneLocal, self.reconstructedSequenceTwoLocal = self.reconstructLocal() 


    def printSimilarityMatrix(self,type="global"):
        if type == "global":
            print('\n'.join(' '.join(map(str, row)) for row in self.matrixGlobal))
        else:
            print('\n'.join(' '.join(map(str, row)) for row in self.matrixLocal))
        return
        
    def printAlignment(self,type="global"):
        if type == "global":
            print(self.reconstructedSequenceOneGlobal + "\n" + self.reconstructedSequenceTwoGlobal)
        else:
            print(self.reconstructedSequenceOneLocal+ "\n" + self.reconstructedSequenceTwoLocal)

    def initMatrixGlobal(self,sequenceOne, sequenceTwo):
        ### height corresponds to sequenceOne, width is sequenceTwo
        ### [ [s2],// [s2],// [s2] ]
        lenSequenceTwo, lenSequenceOne= len(sequenceTwo), len(sequenceOne)
        colLength, rowLength = lenSequenceOne+2, lenSequenceTwo+2
        matrix = [["-"]*rowLength for i in range((colLength))] 
        

        for i in range(lenSequenceOne):
            matrix[i+2][0] = sequenceOne[i]
            matrix[i+2][1] = (i+1)*self.indelPenalty
        for i in range(lenSequenceTwo):
            matrix[0][i+2] = sequenceTwo[i]
            matrix[1][i+2] = (i+1)*self.indelPenalty
        matrix[1][1] = 0
       
        pointers = dict()

        for i in range(1,colLength):
            for j in range(1,rowLength):
                if matrix[i][j] != "-":
                    continue
                
                indelPenalty = max(matrix[i-1][j], matrix[i][j-1]) + self.indelPenalty 
                misMatchPenalty = matrix[i-1][j-1] + self.matchPenalty if matrix[i][0] == matrix[0][j] else matrix[i-1][j-1] + self.mismatchPenalty
                matrix[i][j] = max(indelPenalty, misMatchPenalty)

                if matrix[i][j] == matrix[i-1][j]+ self.indelPenalty :
                    pointers[(i,j)] = [i-1,j]
                elif matrix[i][j] == matrix[i][j-1]+ self.indelPenalty :
                    pointers[(i,j)] = [i,j-1]
                else:
                    pointers[(i,j)] = [i-1,j-1]
        self.globalPointers = pointers
        return matrix

    def reconstructGlobal(self):
        matrix = self.matrixGlobal
        reconstructedOne = ""
        reconstructedTwo = ""
        rowIdx = len(matrix)-1
        colIdx = len(matrix[0])-1
        while True:
            if (rowIdx == 1) and (colIdx == 1):
                break
            rowIdxAbove = rowIdx-1
            colIdxLeft = colIdx-1

            if rowIdxAbove == 0:
                reconstructedTwo = matrix[0][colIdx] + reconstructedTwo
                reconstructedOne = "-" + reconstructedOne
                colIdx -=1
                continue
            if colIdxLeft == 0:
                reconstructedOne = matrix[rowIdx][0] + reconstructedOne
                reconstructedTwo = "-" + reconstructedTwo
                rowIdx -=1
                continue

            newRowIdx, newColIdk = self.globalPointers[(rowIdx,colIdx)]
            
            if newRowIdx != rowIdx and newColIdk == colIdx: 
                reconstructedOne = matrix[rowIdx][0] + reconstructedOne
                reconstructedTwo = "-" + reconstructedTwo
            elif newColIdk != colIdx and newRowIdx == rowIdx:
                reconstructedOne = "-" + reconstructedOne
                reconstructedTwo = matrix[0][colIdx] + reconstructedTwo
            else:
                reconstructedOne = matrix[rowIdx][0] + reconstructedOne
                reconstructedTwo = matrix[0][colIdx] + reconstructedTwo

            rowIdx = newRowIdx
            colIdx = newColIdk

        return reconstructedOne, reconstructedTwo
        
    def initMatrixLocal(self,sequenceOne, sequenceTwo):

        ### height corresponds to sequenceOne, width is sequenceTwo
        ### [ [s2],// [s2],// [s2] ]
        lenSequenceTwo, lenSequenceOne= len(sequenceTwo), len(sequenceOne)
        colLength, rowLength = lenSequenceOne+2, lenSequenceTwo+2
        matrix = [["-"]*rowLength for i in range((colLength))] 

        for i in range(lenSequenceOne):
            matrix[i+2][0] = sequenceOne[i]
            matrix[i+2][1] = 0
        for i in range(lenSequenceTwo):
            matrix[0][i+2] = sequenceTwo[i]
            matrix[1][i+2] = 0
        matrix[1][1] = 0
       
        pointers = dict()
        for i in range(1,colLength):
            for j in range(1,rowLength):
                if matrix[i][j] != "-":
                    continue
                
                indelPenalty = max(matrix[i-1][j], matrix[i][j-1]) + self.indelPenalty 
                misMatchPenalty = matrix[i-1][j-1] + self.matchPenalty if matrix[i][0] == matrix[0][j] else matrix[i-1][j-1] + self.mismatchPenalty
                matrix[i][j] = max(indelPenalty, misMatchPenalty, 0)

                if matrix[i][j] == 0:
                    continue
                elif matrix[i][j] == matrix[i-1][j]+ self.indelPenalty :
                    pointers[(i,j)] = [i-1,j]
                elif matrix[i][j] == matrix[i][j-1]+ self.indelPenalty :
                    pointers[(i,j)] = [i,j-1]
                else:
                    pointers[(i,j)] = [i-1,j-1]
        self.localPointers = pointers

        return matrix

    def reconstructLocal(self):
        matrix = self.matrixLocal
        reconstructedOne=""
        reconstructedTwo=""
        rowLength = len(matrix[0])
        colLength = len(matrix)
        maxEntry = 0
        iMax, jMax = 0 , 0
        for i in range(1,colLength):
            for j in range(1,rowLength):
                if matrix[i][j] > maxEntry:
                    maxEntry = matrix[i][j]
                    iMax, jMax = i, j

        rowIdx,colIdx = iMax, jMax

        while True:
            if matrix[rowIdx][colIdx] == 0:
                break
            rowIdxAbove = rowIdx-1
            colIdxLeft = colIdx-1

            if rowIdxAbove == 0:
                reconstructedTwo = matrix[0][colIdx] + reconstructedTwo
                reconstructedOne = "-" + reconstructedOne
                colIdx -=1
                continue
            if colIdxLeft == 0:
                reconstructedOne = matrix[rowIdx][0] + reconstructedOne
                reconstructedTwo = "-" + reconstructedTwo
                rowIdx -=1
                continue

            newRowIdx, newColIdk = self.localPointers[(rowIdx,colIdx)]
            
            if newRowIdx != rowIdx and newColIdk == colIdx: 
                reconstructedOne = matrix[rowIdx][0] + reconstructedOne
                reconstructedTwo = "-" + reconstructedTwo
            elif newColIdk != colIdx and newRowIdx == rowIdx:
                reconstructedOne = "-" + reconstructedOne
                reconstructedTwo = matrix[0][colIdx] + reconstructedTwo
            else:
                reconstructedOne = matrix[rowIdx][0] + reconstructedOne
                reconstructedTwo = matrix[0][colIdx] + reconstructedTwo

            rowIdx = newRowIdx
            colIdx = newColIdk

        return reconstructedOne, reconstructedTwo
        