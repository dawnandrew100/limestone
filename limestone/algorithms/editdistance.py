from __future__ import annotations
from limestone.algorithms.base import GLOBALBASE as __GLOBALBASE, LOCALBASE as __LOCALBASE
try:
    # external dependencies
    import numpy
    from numpy import float64
    from numpy._typing import NDArray
except ImportError:
    raise ImportError("Please pip install all dependencies from requirements.txt!")

def main():
    qqs = "HOLYWATERISABLESSING"
    sss = "HOLWISBLESSING"

    print(smith_waterman.align(qqs,sss))
    print(smith_waterman.distance(qqs, sss))
    print(smith_waterman.similarity(qqs, sss))
    print(longest_common_subsequence.align(qqs, sss))
    print(longest_common_subsequence.distance(qqs, sss))
    print(longest_common_subsequence.similarity(qqs, sss))

class Wagner_Fischer(__GLOBALBASE): #Levenshtein Distance
    def __init__(self)->None:
        self.gap_penalty = 1

    def __call__(self, querySequence: str, subjectSequence: str)->tuple[NDArray[float64],NDArray[float64]]:
        qs,ss = map(lambda x: x.upper(), [querySequence, subjectSequence])
        qs = [x for x in qs]
        ss = [x for x in ss]
        qs, ss = frontWhiteSpace(qs, ss)

        #matrix initialisation
        self.alignment_score = numpy.zeros((len(qs),len(ss)))
        #pointer matrix to trace optimal alignment
        self.pointer = numpy.zeros((len(qs), len(ss)))
        self.pointer[:,0] = 3
        self.pointer[0,:] = 4
        #initialisation of starter values for first column and first row
        self.alignment_score[:,0] = [n for n in range(len(qs))]
        self.alignment_score[0,:] = [n for n in range(len(ss))]

        for i, query_char in enumerate(qs):
          for j, subject_char in enumerate(ss):
              substitution_cost = 0
              if i == 0 or j == 0:
                  #keeps first row and column consistent throughout all calculations
                  continue
              if query_char != subject_char:
                  substitution_cost = 1
              match = self.alignment_score[i-1][j-1] + substitution_cost
              ugap = self.alignment_score[i-1][j] + self.gap_penalty
              lgap = self.alignment_score[i][j-1] + self.gap_penalty
              tmin = min(match, lgap, ugap)

              self.alignment_score[i][j] = tmin #lowest value is best choice

              if match == tmin: #matrix for traceback based on results from scoring matrix
                  self.pointer[i,j] += 2
              if ugap == tmin:
                  self.pointer[i,j] += 3
              if lgap == tmin:
                  self.pointer[i,j] += 4

        return self.alignment_score, self.pointer

class Lowrance_Wagner(__GLOBALBASE): #Damerau-Levenshtein distance
    def __init__(self)->None:
        self.gap_penalty = 1

    def __call__(self, querySequence: str, subjectSequence: str)->tuple[NDArray[float64],NDArray[float64]]:
        qs,ss = map(lambda x: x.upper(), [querySequence, subjectSequence])
        qs = [x for x in qs]
        ss = [x for x in ss]
        qs, ss = frontWhiteSpace(qs, ss)

        #matrix initialisation
        self.alignment_score = numpy.zeros((len(qs),len(ss)))
        #pointer matrix to trace optimal alignment
        self.pointer = numpy.zeros((len(qs), len(ss)))
        self.pointer[:,0] = 3
        self.pointer[0,:] = 4
        #initialisation of starter values for first column and first row
        self.alignment_score[:,0] = [n for n in range(len(qs))]
        self.alignment_score[0,:] = [n for n in range(len(ss))]

        for i, query_char in enumerate(qs):
          for j, subject_char in enumerate(ss):
              substitution_cost = 0
              if i == 0 or j == 0:
                  #keeps first row and column consistent throughout all calculations
                  continue
              if query_char != subject_char:
                  substitution_cost = 1
              match = self.alignment_score[i-1][j-1] + substitution_cost
              ugap = self.alignment_score[i-1][j] + self.gap_penalty
              lgap = self.alignment_score[i][j-1] + self.gap_penalty
              trans = self.alignment_score[i-2][j-2] + 1 if qs[i] == ss[j-1] and ss[j] == qs[i-1] else float('inf')
              tmin = min(match, lgap, ugap, trans)

              self.alignment_score[i][j] = tmin #lowest value is best choice

              if match == tmin: #matrix for traceback based on results from scoring matrix
                  self.pointer[i,j] += 2
              if ugap == tmin:
                  self.pointer[i,j] += 3
              if lgap == tmin:
                  self.pointer[i,j] += 4
              if trans == tmin:
                  self.pointer[i,j] += 8

        return self.alignment_score, self.pointer

    def align(self, querySequence: str, subjectSequence: str)->str: 
        qs,ss= map(lambda x: x.upper(), [querySequence, subjectSequence])
        _, pointerMatrix = self(qs, ss)

        qs = [x for x in qs]
        ss = [x for x in ss]
        i = len(qs)
        j = len(ss)
        queryAlign= []
        subjectAlign = []

        while i > 0 or j > 0: #looks for match/mismatch/gap starting from bottom right of matrix
          if pointerMatrix[i,j] in [2, 5, 6, 10, 17]:
              #appends match/mismatch then moves to the cell diagonally up and to the left
              queryAlign.append(qs[i-1])
              subjectAlign.append(ss[j-1])
              i -= 1
              j -= 1
          elif pointerMatrix[i,j] in [8, 10, 11, 12, 17]:
              queryAlign.extend([qs[i-1],qs[i-2]])
              subjectAlign.extend([ss[j-2],ss[j-1]])
              i -= 2
              j-= 2
          elif pointerMatrix[i,j] in [3, 5, 7, 11, 17]:
              #appends gap and accompanying nucleotide, then moves to the cell above
              subjectAlign.append('-')
              queryAlign.append(qs[i-1])
              i -= 1
          elif pointerMatrix[i,j] in [4, 6, 7,12, 17]:
              #appends gap and accompanying nucleotide, then moves to the cell to the left
              subjectAlign.append(ss[j-1])
              queryAlign.append('-')
              j -= 1

        queryAlign = "".join(queryAlign[::-1])
        subjectAlign = "".join(subjectAlign[::-1])

        return f"{queryAlign}\n{subjectAlign}"

class Hamming(__GLOBALBASE):
    def __int_pair(self, querySequence: str|int, subjectSequence: str|int) -> bool:
      querySequence, subjectSequence = str(querySequence), str(subjectSequence)
      if querySequence.isalpha() and subjectSequence.isalpha():
          return False
      if querySequence.isdigit() and subjectSequence.isdigit():
          return True 
      raise ValueError("Both sequences must be either all letters or all numbers")
        
    def align(self, querySequence: str|int, subjectSequence: str|int)->str:
        if self.__int_pair(querySequence, subjectSequence):
            qs, ss = int(querySequence), int(subjectSequence)
            return f"{bin(qs)}\n{bin(ss)}"
        return f"{querySequence}\n{subjectSequence}"

    def matrix(self, qs: str, ss: str) -> None:
        return None

    def __call__(self, querySequence: str|int, subjectSequence: str|int)->tuple[int,list[int]]:
        if self.__int_pair(querySequence, subjectSequence):
            qs, ss = bin(querySequence), bin(subjectSequence)
        else:
            qs,ss = map(lambda x: x.upper(), [querySequence,subjectSequence])

        if len(qs) == 1 and len(ss) == 1:
            dist = 1 if qs != ss else 0
            dist_array = [dist]
            return dist, dist_array

        shortlen = min(map(len, [ss,qs]))
        short = qs if len(qs) == shortlen else ss
        long = ss if len(qs) == shortlen else qs

        dist = 0
        dist_array = []

        for i, char in enumerate(short):
            if char != long[i]:
                dist += 1
                dist_array.append(0)
                continue
            dist_array.append(1)

        dist += len(long)-len(short)
        dist_array.extend([1]*(len(long)-len(short)))

        return dist, dist_array

    def distance(self, querySequence: str|int, subjectSequence: str|int)->int:
        if self.__int_pair(querySequence, subjectSequence):
            qs, ss = int(querySequence), int(subjectSequence)
            return bin(qs ^ ss).count("1")
        query = set([(x, y) for (x, y) in enumerate(querySequence)]) 
        subject = set([(x, y) for (x, y) in enumerate(subjectSequence)]) 
        qs,sq = query-subject, subject-query
        dist = max(map(len,[qs,sq]))
        return dist 

    def binary_distance_array(self, querySequence: str, subjectSequence: str)->list[int]:
        _, distarray = self(querySequence, subjectSequence)
        return distarray

    def binary_similarity_array(self, querySequence: str, subjectSequence: str)->list[int]:
        _, distarray = self(querySequence, subjectSequence)
        simarray = [1 if num == 0 else 0 for num in distarray]
        return simarray

class Needleman_Wunsch(__GLOBALBASE):
    def __init__(self, match_score:int = 0, mismatch_penalty:int = 1, gap_penalty:int = 2)->None:
        self.match_score = match_score
        self.mismatch_penalty = mismatch_penalty
        self.gap_penalty = gap_penalty

    def __call__(self, querySequence: str, subjectSequence: str)->tuple[NDArray[float64],NDArray[float64]]:
        qs,ss = map(lambda x: x.upper(), [querySequence, subjectSequence])
        qs = [x for x in qs]
        ss = [x for x in ss]
        qs, ss = frontWhiteSpace(qs, ss)

        #matrix initialisation
        self.alignment_score = numpy.zeros((len(qs),len(ss)))
        #pointer matrix to trace optimal alignment
        self.pointer = numpy.zeros((len(qs), len(ss)))
        self.pointer[:,0] = 3
        self.pointer[0,:] = 4
        #initialisation of starter values for first column and first row
        self.alignment_score[:,0] = [n*self.gap_penalty for n in range(len(qs))]
        self.alignment_score[0,:] = [n*self.gap_penalty for n in range(len(ss))]

        for i, query_char in enumerate(qs):
          for j, subject_char in enumerate(ss):
              if i == 0 or j == 0:
                  #keeps first row and column consistent throughout all calculations
                  continue
              if query_char == subject_char:
                  match = self.alignment_score[i-1][j-1] - self.match_score
              else:
                  match = self.alignment_score[i-1][j-1] + self.mismatch_penalty
              ugap = self.alignment_score[i-1][j] + self.gap_penalty
              lgap = self.alignment_score[i][j-1] + self.gap_penalty
              tmin = min(match, lgap, ugap)

              self.alignment_score[i][j] = tmin #lowest value is best choice

              if match == tmin: #matrix for traceback based on results from scoring matrix
                  self.pointer[i,j] += 2
              if ugap == tmin:
                  self.pointer[i,j] += 3
              if lgap == tmin:
                  self.pointer[i,j] += 4

        return self.alignment_score, self.pointer

class Waterman_Smith_Beyer(__GLOBALBASE):
    def __init__(self, match_score:int = 0, mismatch_penalty:int = 1, new_gap_penalty:int = 1, continue_gap_penalty:int = 1)->None:
        self.match_score = match_score
        self.mismatch_penalty = mismatch_penalty
        self.new_gap_penalty = new_gap_penalty
        self.continue_gap_penalty = continue_gap_penalty

    def __call__(self, querySequence: str, subjectSequence: str)->tuple[NDArray[float64], NDArray[float64]]:
        qs,ss= map(lambda x: x.upper(), [querySequence,subjectSequence])
        qs = [x for x in qs]
        ss = [x for x in ss]
        qs, ss = frontWhiteSpace(qs, ss) 

        #matrix initialisation
        self.alignment_score = numpy.zeros((len(qs),len(ss)))
        #pointer matrix to trace optimal alignment
        self.pointer = numpy.zeros((len(qs), len(ss))) 
        self.pointer[:,0] = 3
        self.pointer[0,:] = 4
        #initialisation of starter values for first column and first row
        self.alignment_score[:,0] = [self.new_gap_penalty + n * self.continue_gap_penalty for n in range(len(qs))]
        self.alignment_score[0,:] = [self.new_gap_penalty + n * self.continue_gap_penalty for n in range(len(ss))] 
        self.alignment_score[0][0] = 0

        for i, subject in enumerate(qs):
            for j, query in enumerate(ss):
                if i == 0 or j == 0:
                    #keeps first row and column consistent throughout all calculations
                    continue
                if subject == query: 
                    matchScore = self.alignment_score[i-1][j-1] - self.match_score
                else:
                    matchScore = self.alignment_score[i-1][j-1] + self.mismatch_penalty
                #both gaps defaulted to continue gap penalty
                ugapScore = self.alignment_score[i-1][j] + self.continue_gap_penalty
                lgapScore = self.alignment_score[i][j-1] + self.continue_gap_penalty
                #if cell before i-1 or j-1 is gap, then this is a gap continuation
                if self.alignment_score[i-1][j] != (self.alignment_score[i-2][j]) + self.new_gap_penalty + self.continue_gap_penalty:
                    ugapScore += self.new_gap_penalty
                if self.alignment_score[i][j-1] != (self.alignment_score[i][j-2]) + self.new_gap_penalty + self.continue_gap_penalty:
                    lgapScore += self.new_gap_penalty
                tmin = min(matchScore, lgapScore, ugapScore)

                self.alignment_score[i][j] = tmin #lowest value is best choice

                #matrix for traceback based on results from scoring matrix
                if matchScore == tmin: 
                    self.pointer[i,j] += 2
                elif ugapScore == tmin:
                    self.pointer[i,j] += 3
                elif lgapScore == tmin:
                    self.pointer[i,j] += 4

        return self.alignment_score, self.pointer

class Hirschberg():
    def __init__(self, match_score: int = 1, mismatch_penalty: int = -1, gap_penalty: int = -2)->None:
        self.match_score = match_score
        self.mismatch_penalty = mismatch_penalty
        self.gap_penalty = gap_penalty

    def __call__(self, querySequence: str, subjectSequence:str)->tuple[str,str]:
        qs,ss = map(lambda x: x.upper(), [querySequence,subjectSequence])

        if len(qs) == 0:
            return '-' * len(ss), ss
        elif len(ss) == 0:
            return qs, '-' * len(qs)
        elif len(qs) == 1 or len(ss) == 1:
            return self._align(qs, ss)
        else:
            xlen = len(qs)
            xmid = xlen // 2
            score_l = self._score(qs[:xmid], ss)
            score_r = self._score(qs[xmid:][::-1], ss[::-1])[::-1]
            ymid = numpy.argmax(score_l + score_r)

            A_left, B_left = self(qs[:xmid], ss[:ymid])
            A_right, B_right = self(qs[xmid:], ss[ymid:])

            return A_left + A_right, B_left + B_right

    def _score(self, qs, ss):
        m, n = len(qs), len(ss)
        score = numpy.zeros((2, n + 1))

        for j in range(1, n + 1):
            score[0][j] = score[0][j - 1] + self.gap_penalty

        for i in range(1, m + 1):
            score[1][0] = score[0][0] + self.gap_penalty
            for j in range(1, n + 1):
                match = score[0][j - 1] + (self.match_score if qs[i - 1] == ss[j - 1] else self.mismatch_penalty)
                delete = score[0][j] + self.gap_penalty
                insert = score[1][j - 1] + self.gap_penalty
                score[1][j] = max(match, delete, insert)
            score[0] = score[1]

        return score[1]

    def _align(self, qs, ss):
        m, n = len(qs), len(ss)
        score = numpy.zeros((m + 1, n + 1))
        pointer = numpy.zeros((m + 1, n + 1))

        for i in range(1, m + 1):
            score[i][0] = score[i - 1][0] + self.gap_penalty
            pointer[i][0] = 1
        for j in range(1, n + 1):
            score[0][j] = score[0][j - 1] + self.gap_penalty
            pointer[0][j] = 2

        for i in range(1, m + 1):
            for j in range(1, n + 1):
                if qs[i-1] == ss[j-1]: 
                    match = score[i-1][j-1] + self.match_score
                else:
                    match = score[i-1][j-1] + self.mismatch_penalty
                delete = score[i - 1][j] + self.gap_penalty
                insert = score[i][j - 1] + self.gap_penalty
                score[i][j] = max(match, delete, insert)
                if score[i][j] == match:
                    pointer[i][j] = 3
                elif score[i][j] == delete:
                    pointer[i][j] = 1
                else:
                    pointer[i][j] = 2

        queryAlign, subjectAlign = "", ""
        i, j = m, n

        while i > 0 or j > 0:
            if pointer[i][j] == 3:
                queryAlign = qs[i - 1] + queryAlign
                subjectAlign = ss[j - 1] + subjectAlign
                i -= 1
                j -= 1
            elif pointer[i][j] == 1:
                queryAlign = qs[i - 1] + queryAlign
                subjectAlign = '-' + subjectAlign
                i -= 1
            else:
                queryAlign = '-' + queryAlign
                subjectAlign = ss[j - 1] + subjectAlign
                j -= 1

        return queryAlign, subjectAlign

    def align(self, qs: str, ss: str) -> str:
        queryAlign, subjectAlign = self(qs, ss)
        return f"{queryAlign}\n{subjectAlign}"

class Smith_Waterman(__LOCALBASE):
    def __init__(self, match_score:int = 1, mismatch_penalty:int = 1, gap_penalty:int = 2)->None:
        self.match_score = match_score
        self.mismatch_penalty = mismatch_penalty
        self.gap_penalty = gap_penalty

    def __call__(self, querySequence: str, subjectSequence: str)-> NDArray[float64]: 
        qs,ss = map(lambda x: x.upper(), [querySequence,subjectSequence])
        qs = [x for x in qs]
        ss = [x for x in ss]
        qs, ss = frontWhiteSpace(qs, ss) 

        #matrix initialisation
        self.alignment_score = numpy.zeros((len(qs),len(ss))) 

        for i, query_char in enumerate(qs):
          for j, subject_char in enumerate(ss):
            if j == 0 or i == 0:
                #keeps first row and column consistent throughout all calculations
                continue

            if query_char == subject_char: 
                match = self.alignment_score[i-1][j-1] + self.match_score
            else:
                match = self.alignment_score[i-1][j-1] - self.mismatch_penalty
 
            ugap = self.alignment_score[i-1][j] - self.gap_penalty 
            lgap = self.alignment_score[i][j-1] - self.gap_penalty 
            tmax = max(0, match, lgap, ugap) 
            self.alignment_score[i][j] = tmax

        return self.alignment_score
 
    def align(self, querySequence: str, subjectSequence: str)->str:
      qs,ss = map(lambda x: x.upper(), [querySequence, subjectSequence])
      matrix = self(qs, ss)

      qs = [x for x in qs]
      ss = [x for x in ss]

      if matrix.max() == 0:
        return "There is no local alignment!"

      #finds the largest value closest to bottom right of matrix
      i, j = list(numpy.where(matrix == matrix.max()))
      i, j = i[-1], j[-1]

      subjectAlign = []
      queryAlign= []
      score = matrix.max()

      while score > 0:
          score = matrix[i][j]
          if score == 0:
              break
          queryAlign.append(qs[i-1])
          subjectAlign.append(ss[j-1])
          i -= 1
          j -= 1

      queryAlign = "".join(queryAlign[::-1])
      subjectAlign = "".join(subjectAlign[::-1])

      return f"{queryAlign}\n{subjectAlign}"

class Longest_Common_Subsequence(__LOCALBASE):
    def __init__(self):
        self.match_score = 1

    def __call__(self, querySequence: str, subjectSequence: str)-> NDArray[float64]: 
        qs,ss = map(lambda x: x.upper(), [querySequence,subjectSequence])
        qs = [x for x in qs]
        ss = [x for x in ss]
        qs, ss = frontWhiteSpace(qs, ss) 

        #matrix initialisation
        self.alignment_score = numpy.zeros((len(qs),len(ss))) 

        for i, query_char in enumerate(qs):
          for j, subject_char in enumerate(ss):
            if j == 0 or i == 0:
                #keeps first row and column consistent throughout all calculations
                continue

            if query_char == subject_char: 
                match = self.alignment_score[i-1][j-1] + self.match_score
            else:
                match = max(self.alignment_score[i][j-1], self.alignment_score[i-1][j]) 
 
            self.alignment_score[i][j] = match

        return self.alignment_score
 
    def align(self, querySequence: str, subjectSequence: str)->str:
      qs,ss = map(lambda x: x.upper(), [querySequence, subjectSequence])
      matrix = self(qs, ss)

      qs = [x for x in qs]
      ss = [x for x in ss]

      if matrix.max() == 0:
        return "There is no common subsequence!"

      #starts backtrack in bottom right of matrix
      i, j = len(querySequence), len(subjectSequence)

      subjectAlign = []
      queryAlign= []

      while matrix[i, j] > 0:
          if i == 0 and j == 0:
              break
          if querySequence[i-1] == subjectSequence[j-1]:
              queryAlign.append(qs[i-1])
              subjectAlign.append(ss[j-1])
              i -= 1
              j -= 1
          elif matrix[i-1,j] >= matrix[i, j-1]:
              i -= 1
          elif matrix[i, j-1] >= matrix[i-1, j]:
              j -= 1

      queryAlign = "".join(queryAlign[::-1])
      subjectAlign = "".join(subjectAlign[::-1])

      return f"{queryAlign}\n{subjectAlign}"



def frontWhiteSpace(qs: list[str], ss: list[str])->tuple[list[str],list[str]]: 
    #adds leading white space so that matrix works
    qs = rjustlist(qs)
    ss = rjustlist(ss)
    return qs, ss

def ljustlist(sequence: list[str], n: int, fillvalue='')->list[str]:
    return sequence + [fillvalue] * (n - len(sequence)) 

def rjustlist(sequence: list[str], fillvalue='')->list[str]:
    return [fillvalue] + sequence


hamming = Hamming()
wagner_fischer = Wagner_Fischer()
needleman_wunsch = Needleman_Wunsch()
waterman_smith_beyer = Waterman_Smith_Beyer()
smith_waterman = Smith_Waterman()
hirschberg = Hirschberg()
lowrance_wagner = Lowrance_Wagner()
longest_common_subsequence = Longest_Common_Subsequence()

if __name__ == "__main__":
    main()
