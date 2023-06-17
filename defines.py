

class SequenceAlignment:
    def __init__(self,sequenceOne,sequenceTwo,type):
        if type == "global":
            self.matrix = self.initMatrixGlobal(sequenceOne,sequenceTwo)

    def __str__(self):
        return '\n'.join(' '.join(map(str, row)) for row in self.matrix)


    def initMatrixGlobal(self,sequenceOne, sequenceTwo):
        ### height corresponds to sequenceOne, width is sequenceTwo
        ### [ [s2],// [s2],// [s2] ]
        lenSequenceTwo, lenSequenceOne= len(sequenceTwo), len(sequenceOne)
        maxLength,minLength = max(lenSequenceOne,lenSequenceTwo), min(lenSequenceOne,lenSequenceTwo)
        colLength, rowLength = lenSequenceOne+2, lenSequenceTwo+2
        matrix = [["-"]*rowLength for i in range((colLength))] 
        

        for i in range(lenSequenceOne):
            matrix[i+2][0] = sequenceOne[i]
            matrix[i+2][1] = -(i+1)
        for i in range(lenSequenceTwo):
            matrix[0][i+2] = sequenceTwo[i]
            matrix[1][i+2] = -(i+1)
        matrix[1][1] = 0
       

        for i in range(1,colLength):
            for j in range(1,rowLength):
                if matrix[i][j] != "-":
                    continue
                
                indelPenalty = max(matrix[i-1][j], matrix[i][j-1])-2
                misMatchPenalty = matrix[i-1][j-1] + 1 if matrix[i][0] == matrix[0][j] else matrix[i-1][j-1]-1
                matrix[i][j] = max(indelPenalty, misMatchPenalty)

        return matrix

    def reconstructGlobal(self):
        matrix = self.matrix
        reconstructedOne=""
        reconstructedTwo=""
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

            maxChar = max(matrix[rowIdxAbove][colIdx],matrix[rowIdx][colIdxLeft], matrix[rowIdxAbove][colIdxLeft])
            
            if matrix[rowIdxAbove][colIdxLeft] == maxChar:
                reconstructedOne = matrix[rowIdx][0] + reconstructedOne
                reconstructedTwo = matrix[0][colIdx] + reconstructedTwo
                colIdx -= 1
                rowIdx -= 1
            elif matrix[rowIdx][colIdxLeft] == maxChar:
                reconstructedOne = "-" + reconstructedOne
                reconstructedTwo = matrix[0][colIdx] + reconstructedTwo
                colIdx -= 1
            else:
                reconstructedOne = matrix[rowIdx][0] + reconstructedOne
                reconstructedTwo = "-" + reconstructedTwo
                rowIdx -= 1
            



        return reconstructedOne, reconstructedTwo


        

