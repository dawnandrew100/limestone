from __future__ import annotations
from re import sub
try:
    # external dependency
    import numpy
    from numpy import float64
    from numpy._typing import NDArray
except ImportError:
    numpy = None  

def main():
    qqs = "AT"
    sss = "AAGT"

    print(waterman_smith_beyer.matrix(qqs, sss))

class _GLOBALBASE():
  def matrix(self, querySequence: str, subjectSequence: str)->list[list[float]]:
    matrix, _ = self(querySequence, subjectSequence)
    return matrix

  def distance(self, querySequence: str, subjectSequence: str)->float:
    matrix, _ = self(querySequence, subjectSequence)
    return matrix[matrix.shape[0]-1,matrix.shape[1]-1]

  def similarity(self, querySequence: str, subjectSequence: str)->float:
    return max(len(querySequence),len(subjectSequence)) - self.distance(querySequence, subjectSequence)

  def normalized_distance(self, querySequence: str, subjectSequence: str)->float:
    dist = self.distance(querySequence, subjectSequence)
    return dist/max(map(len, [querySequence,subjectSequence]))

  def normalized_similarity(self, querySequence: str, subjectSequence: str)->float:
    sim = self.similarity(querySequence, subjectSequence)
    return sim/max(map(len, [querySequence,subjectSequence]))

  def align(self, querySequence: str, subjectSequence: str)->str: 
      qs,ss= map(lambda x: x.upper(), [querySequence,subjectSequence])
      _, pointerMatrix = self(qs, ss)
      
      qs = [x for x in qs]
      ss = [x for x in ss]
      i = len(qs)
      j = len(ss)
      queryAlign= []
      subjectAlign = []

      while i > 0 or j > 0: #looks for match/mismatch/gap starting from bottom right of matrix
        if pointerMatrix[i,j] in [2, 5, 6, 9]:
            #appends match/mismatch then moves to the cell diagonally up and to the left
            queryAlign.append(qs[i-1])
            subjectAlign.append(ss[j-1])
            i -= 1
            j -= 1
        elif pointerMatrix[i,j] in [3, 5, 7, 9]:
          #appends gap and accompanying nucleotide, then moves to the cell above
            subjectAlign.append('-')
            queryAlign.append(qs[i-1])
            i -= 1
        elif pointerMatrix[i,j] in [4, 6, 7, 9]:
            #appends gap and accompanying nucleotide, then moves to the cell to the left
            subjectAlign.append(ss[j-1])
            queryAlign.append('-')
            j -= 1

      queryAlign = "".join(queryAlign[::-1])
      subjectAlign = "".join(subjectAlign[::-1])

      return f"{queryAlign}\n{subjectAlign}"

class _LOCALBASE():
  #All local base functions currently only used by Smith Waterman
  def matrix(self, querySequence: str, subjectSequence: str)->list[list[float]]:
    matrix = self(querySequence, subjectSequence)
    return matrix

  def similarity(self, querySequence: str,subjectSequence: str)->float:
    matrix  = self(querySequence, subjectSequence)
    return matrix.max()

  def distance(self, querySequence: str, subjectSequence: str)->float:
    return max(map(len, [querySequence,subjectSequence])) - self.similarity(querySequence, subjectSequence)

  def normalized_distance(self, querySequence: str, subjectSequence: str)->float:
    dist = self.distance(querySequence, subjectSequence)
    return dist/max(map(len, [querySequence,subjectSequence]))

  def normalized_similarity(self, querySequence: str, subjectSequence: str)->float:
    similarity = self.similarity(querySequence, subjectSequence)
    return similarity/max(map(len, [querySequence,subjectSequence]))

  def align(self, querySequence: str, subjectSequence: str)->str:
    qs,ss = map(lambda x: x.upper(), [querySequence,subjectSequence])
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
        queryAlign.append(qs[j-1])
        subjectAlign.append(ss[i-1])
        i -= 1
        j -= 1

    queryAlign = "".join(queryAlign[::-1])
    subjectAlign = "".join(subjectAlign[::-1])

    return f"{queryAlign}\n{subjectAlign}"

class Hamming(_GLOBALBASE):
    def align(self, querySequence: str, subjectSequence: str)->str:
        return f"{querySequence}\n{subjectSequence}"

    def matrix(self, qs: str, ss: str) -> None:
        return None

    def __call__(self, querySequence: str, subjectSequence: str)->tuple[int,list[int]]:
      if not numpy:
          raise ImportError('Please pip install numpy!')
      
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

    def distance(self, querySequence: str, subjectSequence: str)->int:
        query = set([(x,y) for (x,y) in enumerate(querySequence)]) 
        subject = set([(x,y) for (x,y) in enumerate(subjectSequence)]) 
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

class Needleman_Wunsch(_GLOBALBASE):
  def __init__(self, match_score:int = 0, mismatch_penalty:int = 1, gap_penalty:int = 2)->None:
    self.match_score = match_score
    self.mismatch_penalty = mismatch_penalty
    self.gap_penalty = gap_penalty

  def __call__(self, querySequence: str, subjectSequence: str)->tuple[NDArray[float64],NDArray[float64]]:
      if not numpy:
          raise ImportError('Please pip install numpy!')

      qs,ss = map(lambda x: x.upper(), [querySequence,subjectSequence])
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

class Levenshtein(Needleman_Wunsch):
  def __init__(self):
    self.match_score = 0
    self.mismatch_penalty = 1
    self.gap_penalty = 1

class Smith_Waterman(_LOCALBASE):
  def __init__(self, match_score:int = 1, mismatch_penalty:int = 1, gap_penalty:int = 2)->None:
    self.match_score = match_score
    self.mismatch_penalty = mismatch_penalty
    self.gap_penalty = gap_penalty

  def __call__(self, subjectSequence: str, querySequence: str)-> NDArray[float64]: 
      if not numpy:
          raise ImportError('Please pip install numpy!')

      qs,ss = map(lambda x: x.upper(), [querySequence,subjectSequence])
      qs = [x for x in qs]
      ss = [x for x in ss]
      qs, ss = frontWhiteSpace(qs, ss) 
      
      #matrix initialisation
      self.alignment_score = numpy.zeros((len(qs),len(ss))) 

      self.alignment_score[:,0] = 0
      self.alignment_score[0,:] = 0
      
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

class Waterman_Smith_Beyer(_GLOBALBASE):
  def __init__(self, match_score:int = 0, mismatch_penalty:int = 1, new_gap_penalty:int = 3, continue_gap_penalty:int = 1)->None:
      self.match_score = match_score
      self.mismatch_penalty = mismatch_penalty
      self.new_gap_penalty = new_gap_penalty
      self.continue_gap_penalty = continue_gap_penalty

  def __call__(self, querySequence: str, subjectSequence: str)->tuple[NDArray[float64], NDArray[float64]]:
      if not numpy:
          raise ImportError('Please pip install numpy!')

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
        if not numpy:
            raise ImportError('Please pip install numpy!')

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

def frontWhiteSpace(qs: list[str], ss: list[str])->tuple[list[str],list[str]]: 
    #adds leading white space so that matrix works
    qs = rjustlist(qs)
    ss = rjustlist(ss)
    return qs, ss

def ljustlist(sequence: list[str], n: int, fillvalue='')->list[str]:
    return sequence + [fillvalue] * (n - len(sequence)) 

def rjustlist(sequence: list[str], fillvalue='')->list[str]:
    return [fillvalue] + sequence


needleman_wunsch = Needleman_Wunsch()
waterman_smith_beyer = Waterman_Smith_Beyer()
smith_waterman = Smith_Waterman()
levenshtein = Levenshtein()
hirschberg = Hirschberg()
hamming = Hamming()

if __name__ == "__main__":
    main()
