from __future__ import annotations
try:
    # external dependency
    import numpy
except ImportError:
    numpy = None  

def main():
    print(hirschbergDist.align("AGTACGCA","TATGC"))
    print(needlemanWunsch.align("AGTACGCA","TATGC"))

def ljustlist(sequence: list[str], n: int, fillvalue='')->list[str]:
  return sequence + [fillvalue] * (n - len(sequence)) 

def rjustlist(sequence: list[str], fillvalue='')->list[str]:
  return [fillvalue] + sequence

class _BASE():
  @staticmethod
  def maximum(querySequence: str|list[str], subjectSequence: str|list[str])->int:
      sequences = [querySequence,subjectSequence]
      return max(map(len, sequences))

class _GLOBALBASE(_BASE):  
  def scoreMatrix(self, querySequence: str|list[str], subjectSequence: str|list[str])->list[list[int]]:
    matrix, _ = self(querySequence, subjectSequence)
    return matrix
  
  def distance(self, querySequence: str|list[str], subjectSequence: str|list[str])->int:
    matrix, _ = self(querySequence, subjectSequence)
    return matrix[matrix.shape[0]-1,matrix.shape[1]-1]

  def similarity(self, querySequence: str|list[str], subjectSequence: str|list[str])->int:
    return max(len(subjectSequence),len(querySequence)) - self.distance(querySequence, subjectSequence)

  def normalized_distance(self, querySequence: str|list[str], subjectSequence: str|list[str])->int:
    dist = self.distance(querySequence, subjectSequence)
    sequences = [querySequence,subjectSequence]
    return dist/max(map(len, sequences))

  def normalized_similarity(self, querySequence: str|list[str], subjectSequence: str|list[str])->int:
    sim = self.similarity(querySequence, subjectSequence)
    sequences = [querySequence,subjectSequence]
    return sim/max(map(len, sequences))

  def align(self, querySequence: str|list[str], subjectSequence: str|list[str])->str: 
      #currently only used for Needleman Wunsch
      if isinstance(subjectSequence, str):
          querySequence,subjectSequence = map(lambda x: x.upper(), [querySequence,subjectSequence])
          querySequence = [x for x in querySequence]
          subjectSequence = [x for x in subjectSequence]
      
      _, pointerMatrix = self(querySequence, subjectSequence)
  
      subjectSequence, querySequence = frontWhiteSpace(subjectSequence, querySequence)
      i = len(querySequence)+1
      j = len(subjectSequence)+1
      subjectAlign = ""
      queryAlign= ""
      subjectSequence, querySequence = frontWhiteSpace(subjectSequence, querySequence)
      
      while i > 0 or j > 0: #looks for match/mismatch/gap starting from bottom right of matrix
        if pointerMatrix[i,j] in [2, 5, 6, 9]:
            #appends match/mismatch then moves to the cell diagonally up and to the left
            queryAlign = querySequence[j-1] + queryAlign
            subjectAlign = subjectSequence[i-1] + subjectAlign
            i -= 1
            j -= 1
        elif pointerMatrix[i,j] in [3, 5, 7, 9]:
            #appends gap and accompanying nucleotide, then moves to the cell to the left
            subjectAlign = subjectSequence[i-1] +subjectAlign
            queryAlign = '-' +queryAlign
            i -= 1
        elif pointerMatrix[i,j] in [4, 6, 7, 9]:
            #appends gap and accompanying nucleotide, then moves to the cell above
            subjectAlign = '-' + subjectAlign
            queryAlign = querySequence[j-1] + queryAlign
            j -= 1
      return f"{queryAlign}\n{subjectAlign}"

class _LOCALBASE(_BASE):
  #All local base functions currently only used by Smith Waterman
  def scoreMatrix(self, querySequence: str|list[str], subjectSequence: str|list[str])->list[list[int]]:
    matrix = self(querySequence, subjectSequence)
    return matrix
  
  def similarity(self, querySequence: str|list[str],subjectSequence: str|list[str])->int:
    scoreMatrix  = self(querySequence, subjectSequence)
    return scoreMatrix.max()

  def distance(self, querySequence: str|list[str], subjectSequence: str|list[str])->int:
    sequences = [querySequence,subjectSequence]
    return max(map(len, sequences)) - self.similarity(querySequence, subjectSequence)

  def normalized_distance(self, querySequence: str|list[str], subjectSequence: str|list[str])->int:
    dist = self.distance(querySequence, subjectSequence)
    sequences = [querySequence,subjectSequence]
    return dist/max(map(len, sequences))

  def normalized_similarity(self, querySequence: str|list[str], subjectSequence: str|list[str])->int:
    similarity = self.similarity(querySequence, subjectSequence)
    sequences = [querySequence,subjectSequence]
    return similarity/max(map(len, sequences))

  def align(self, querySequence: str|list[str], subjectSequence: str|list[str])->str:
    querySequence,subjectSequence = map(lambda x: x.upper(), [querySequence,subjectSequence])
    scoreMatrix = self(querySequence, subjectSequence)

    querySequence = [x for x in querySequence]
    subjectSequence = [x for x in subjectSequence]

    if scoreMatrix.max() == 0:
      return "There is no local alignment!"

    #finds the largest value closest to bottom right of matrix
    i, j = list(numpy.where(scoreMatrix == scoreMatrix.max()))
    i, j = i[-1], j[-1]

    subjectAlign = ""
    queryAlign= ""
    score = scoreMatrix.max()

    while score > 0:
        score = scoreMatrix[i][j]
        if score == 0:
            break
        queryAlign = querySequence[j-1] + queryAlign
        subjectAlign = subjectSequence[i-1] + subjectAlign
        i -= 1
        j -= 1

    return f"{queryAlign}\n{subjectAlign}"
    
class needleman_wunsch(_GLOBALBASE):
  def __init__(self, match_score:int = 0, mismatch_penalty:int = 1, gap_penalty:int = 2)->None:
    self.match_score = match_score
    self.mismatch_penalty = mismatch_penalty
    self.gap_penalty = gap_penalty

  def __call__(self, querySequence: str|list[str], subjectSequence: str|list[str])->tuple[int,int]:
      if not numpy:
          raise ImportError('Please pip install numpy!')
      
      if isinstance(subjectSequence, str):
          querySequence,subjectSequence = map(lambda x: x.upper(), [querySequence,subjectSequence])
          querySequence = [x for x in querySequence]
          subjectSequence = [x for x in subjectSequence]

      #This needs to be called twice to work (will break if done any other way)
      subjectSequence, querySequence = frontWhiteSpace(subjectSequence, querySequence) 
      subjectSequence, querySequence = frontWhiteSpace(subjectSequence, querySequence) 

      #matrix initialisation
      self.alignment_score = numpy.zeros((len(subjectSequence),len(querySequence)))
      #pointer matrix to trace optimal alignment
      self.pointer = numpy.zeros((len(subjectSequence)+1, len(querySequence)+1)) 
      self.pointer[:,0] = 3
      self.pointer[0,:] = 4
      #first row and column of matrix consisting of gap scores
      init1 = [n*self.gap_penalty for n in range(len(subjectSequence))]
      init2 = [n*self.gap_penalty for n in range(len(querySequence))]

      for i, _ in enumerate(subjectSequence):
        for j, _ in enumerate(querySequence):
          #keeps first row and column consistent throughout all calculations
          self.alignment_score[:,0] = init1
          self.alignment_score[:1,] = init2

          if subjectSequence[i] == querySequence[j]: 
              match = self.alignment_score[i-1][j-1] - self.match_score
          else:
              match = self.alignment_score[i-1][j-1] + self.mismatch_penalty

          lgap = self.alignment_score[i][j-1] + self.gap_penalty 
          ugap = self.alignment_score[i-1][j] + self.gap_penalty 
          tmin = min(match, lgap, ugap)

          self.alignment_score[i][j] = tmin #lowest value is best choice

          if match == tmin: #matrix for traceback based on results from scoring matrix
            self.pointer[i+1,j+1] += 2
          if ugap == tmin:
            self.pointer[i+1,j+1] += 3
          if lgap == tmin:
            self.pointer[i+1,j+1] += 4

      return self.alignment_score, self.pointer
    
class levenshtein(needleman_wunsch):
  def __init__(self):
    self.match_score = 0
    self.mismatch_penalty = 1
    self.gap_penalty = 1
      
class smith_waterman(_LOCALBASE):
  def __init__(self, match_score:int = 1, mismatch_penalty:int = 1, gap_penalty:int = 2)->None:
    self.match_score = match_score
    self.mismatch_penalty = mismatch_penalty
    self.gap_penalty = gap_penalty

  def __call__(self, subjectSequence: str|list[str], querySequence: str|list[str])-> int: 
      if not numpy:
          raise ImportError('Please pip install numpy!')
      
      if isinstance(subjectSequence, str):
          querySequence,subjectSequence = map(lambda x: x.upper(), [querySequence,subjectSequence])
          querySequence = [x for x in querySequence]
          subjectSequence = [x for x in subjectSequence]

      subjectSequence, querySequence = frontWhiteSpace(subjectSequence, querySequence) 
      #matrix initialisation
      self.alignment_score = numpy.zeros((len(subjectSequence),len(querySequence))) 

      for i, _ in enumerate(subjectSequence):
        for j, _ in enumerate(querySequence):
          #keeps first row and column consistent throughout all calculations
          self.alignment_score[:,0] = 0
          self.alignment_score[:1,] = 0

          if subjectSequence[i] == querySequence[j]: 
              match = self.alignment_score[i-1][j-1] + self.match_score
          else:
              match = self.alignment_score[i-1][j-1] - self.mismatch_penalty
  
          lgap = self.alignment_score[i][j-1] - self.gap_penalty 
          ugap = self.alignment_score[i-1][j] - self.gap_penalty 
          tmax = max(0, match, lgap, ugap) 
  
          self.alignment_score[i][j] = tmax
    
      return self.alignment_score

class waterman_smith_beyer(_GLOBALBASE):
  def __init__(self, match_score:int = 0, mismatch_penalty:int = 1, new_gap_penalty:int = 3, continue_gap_penalty:int = 1)->None:
      self.match_score = match_score
      self.mismatch_penalty = mismatch_penalty
      self.new_gap_penalty = new_gap_penalty
      self.continue_gap_penalty = continue_gap_penalty

  def __call__(self, querySequence: str|list[str], subjectSequence: str|list[str])->tuple[int, int]:
      if not numpy:
          raise ImportError('Please pip install numpy!')
      
      if isinstance(subjectSequence, str):
          querySequence,subjectSequence = map(lambda x: x.upper(), [querySequence,subjectSequence])
          querySequence = [x for x in querySequence]
          subjectSequence = [x for x in subjectSequence]
      
      #This needs to be called twice to work (will break if done any other way)
      subjectSequence, querySequence = frontWhiteSpace(subjectSequence, querySequence) 
      subjectSequence, querySequence = frontWhiteSpace(subjectSequence, querySequence)
      
      #matrix initialisation
      self.alignment_score = numpy.zeros((len(subjectSequence),len(querySequence)))
      #pointer matrix to trace optimal alignment
      self.pointer = numpy.zeros((len(subjectSequence)+1, len(querySequence)+1)) 
      self.pointer[:,0] = 3
      self.pointer[0,:] = 4
      #first row and column of matrix consisting of gap scores
      init1 = [n*self.continue_gap_penalty+self.new_gap_penalty for n in range(len(subjectSequence))]
      init2 = [n*self.continue_gap_penalty+self.new_gap_penalty for n in range(len(querySequence))]
  
      for i, _ in enumerate(subjectSequence):
        for j, _ in enumerate(querySequence):
          #keeps first row and column consistent throughout all calculations
          self.alignment_score[:,0] = init1
          self.alignment_score[:1,] = init2
          self.alignment_score[0][0] = 0
  
          if subjectSequence[i] == querySequence[j]: 
            matchScore = self.alignment_score[i-1][j-1] - self.match_score
          else:
            matchScore = self.alignment_score[i-1][j-1] + self.mismatch_penalty
          #both gaps defaulted to new gap penalty
          lgapScore = self.alignment_score[i][j-1] + self.new_gap_penalty + self.continue_gap_penalty
          ugapScore = self.alignment_score[i-1][j] + self.new_gap_penalty + self.continue_gap_penalty
          #if cell before i-1 or j-1 is gap, then this is a gap continuation
          if self.alignment_score[i][j-1] == (self.alignment_score[i][j-2]) + self.new_gap_penalty + self.continue_gap_penalty:
            lgapScore = lgapScore - self.new_gap_penalty
          if self.alignment_score[i-1][j] == (self.alignment_score[i-2][j]) + self.new_gap_penalty + self.continue_gap_penalty:
            ugapScore = ugapScore - self.new_gap_penalty
          tmin = min(matchScore, lgapScore, ugapScore)
  
          self.alignment_score[i][j] = tmin #lowest value is best choice
  
          if matchScore == tmin: #matrix for traceback based on results from scoring matrix
            self.pointer[i+1,j+1] += 2
          if ugapScore == tmin:
            self.pointer[i+1,j+1] += 3
          if lgapScore == tmin:
            self.pointer[i+1,j+1] += 4
  
      return self.alignment_score, self.pointer

class hirschberg(_GLOBALBASE):
    def __init__(self, match_score: int = 1, mismatch_penalty: int = -1, gap_penalty: int = -2):
        self.match_score = match_score
        self.mismatch_penalty = mismatch_penalty
        self.gap_penalty = gap_penalty
    
    def __call__(self, querySequence, subjectSequence):
        if isinstance(querySequence, str):
            querySequence = querySequence.upper()
            subjectSequence = subjectSequence.upper()
        
        if len(querySequence) == 0:
            return '-' * len(subjectSequence), subjectSequence
        elif len(subjectSequence) == 0:
            return querySequence, '-' * len(querySequence)
        elif len(querySequence) == 1 or len(subjectSequence) == 1:
            return self._nw_align(querySequence, subjectSequence)
        else:
            xlen = len(querySequence)
            xmid = xlen // 2
            score_l = self._nw_score(querySequence[:xmid], subjectSequence)
            score_r = self._nw_score(querySequence[xmid:][::-1], subjectSequence[::-1])[::-1]
            ymid = numpy.argmax(score_l + score_r)

            A_left, B_left = self(querySequence[:xmid], subjectSequence[:ymid])
            A_right, B_right = self(querySequence[xmid:], subjectSequence[ymid:])
            return A_left + A_right, B_left + B_right

    def _nw_score(self, querySequence, subjectSequence):
        m, n = len(querySequence), len(subjectSequence)
        score = numpy.zeros((2, n + 1))

        for j in range(1, n + 1):
            score[0][j] = score[0][j - 1] + self.gap_penalty

        for i in range(1, m + 1):
            score[1][0] = score[0][0] + self.gap_penalty
            for j in range(1, n + 1):
                match = score[0][j - 1] + (self.match_score if querySequence[i - 1] == subjectSequence[j - 1] else self.mismatch_penalty)
                delete = score[0][j] + self.gap_penalty
                insert = score[1][j - 1] + self.gap_penalty
                score[1][j] = max(match, delete, insert)
            score[0] = score[1]

        return score[1]

    def _nw_align(self, querySequence, subjectSequence):
        m, n = len(querySequence), len(subjectSequence)
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
                match = score[i - 1][j - 1] + (self.match_score if querySequence[i - 1] == subjectSequence[j - 1] else self.mismatch_penalty)
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
                queryAlign = querySequence[i - 1] + queryAlign
                subjectAlign = subjectSequence[j - 1] + subjectAlign
                i -= 1
                j -= 1
            elif pointer[i][j] == 1:
                queryAlign = querySequence[i - 1] + queryAlign
                subjectAlign = '-' + subjectAlign
                i -= 1
            else:
                queryAlign = '-' + queryAlign
                subjectAlign = subjectSequence[j - 1] + subjectAlign
                j -= 1

        return queryAlign, subjectAlign

    def align(self, querySequence: str | list[str], subjectSequence: str | list[str]) -> str:
        aligned_query, aligned_subject = self(querySequence, subjectSequence)
        return f"{aligned_query}\n{aligned_subject}"

def frontWhiteSpace(querySequence: list[str], subjectSequence: list[str])->tuple[list[str],list[str]]: 
    #adds leading white space so that matrix works
    subjectSequence = rjustlist(subjectSequence)
    querySequence = rjustlist(querySequence)

    return subjectSequence, querySequence

needlemanWunsch = needleman_wunsch()
watermanSmithBeyer = waterman_smith_beyer()
smithWaterman = smith_waterman()
levenshteinDist = levenshtein()
hirschbergDist = hirschberg()

if __name__ == "__main__":
    main()
