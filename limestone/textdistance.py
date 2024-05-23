from __future__ import annotations
try:
    # external dependency
    import numpy
except ImportError:
    numpy = None  

def ljustlist(sequence: list[str], n: int, fillvalue='')->list[str]:
  return sequence + [fillvalue] * (n - len(sequence)) 

def rjustlist(sequence: list[str], fillvalue='')->list[str]:
  return [fillvalue] + sequence

class _BASE():
  @staticmethod
  def maximum(querySequence: str|list[str], subjectSequence: str|list[str])->int:
      sequences = [querySequence,subjectSequence]
      return max(map(len, sequences))

  def scoreMatrix(self, querySequence: str|list[str], subjectSequence: str|list[str])->list[list[int]]:
    matrix, _ = self(querySequence, subjectSequence)
    return matrix

class _GLOBALBASE(_BASE):  
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
      subjectAlign = []
      queryAlign= []
      subjectSequence, querySequence = frontWhiteSpace(subjectSequence, querySequence)
      
      while i > 0 or j > 0: #looks for match/mismatch/gap starting from bottom right of matrix
        if pointerMatrix[i,j] in [2, 5, 6, 9]:
            #appends match/mismatch then moves to the cell diagonally up and to the left
            queryAlign.append(querySequence[j-1])
            subjectAlign.append(subjectSequence[i-1])
            i -= 1
            j -= 1
        elif pointerMatrix[i,j] in [3, 5, 7, 9]:
            #appends gap and accompanying nucleotide, then moves to the cell to the left
            subjectAlign.append(subjectSequence[i-1])
            queryAlign.append('-')
            i -= 1
        elif pointerMatrix[i,j] in [4, 6, 7, 9]:
            #appends gap and accompanying nucleotide, then moves to the cell above
            subjectAlign.append('-')
            queryAlign.append(querySequence[j-1])
            j -= 1
      # Reverses the strings
      subjectAlign = ''.join(subjectAlign)[::-1]
      queryAlign= ''.join(queryAlign)[::-1]
      '\n'.join([subjectAlign, queryAlign])
      return f"{queryAlign}\n{subjectAlign}"

class _LOCALBASE(_BASE):
  #All local base functions currently only used by Smith Waterman
  def similarity(self, querySequence: str|list[str],subjectSequence: str|list[str])->int:
    scoreMatrix, _ = self(querySequence, subjectSequence)
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
    scoreMatrix, pointerMatrix = self(querySequence, subjectSequence)

    querySequence = [x for x in querySequence]
    subjectSequence = [x for x in subjectSequence]

    if scoreMatrix.max() == 0:
      return "There is no local alignment!"

    #finds the largest value closest to bottom right of matrix
    i, j = list(numpy.where(scoreMatrix == scoreMatrix.max()))
    i, j = i[-1], j[-1]

    subjectAlign = []
    queryAlign= []
    qtemp = ''
    stemp = ''

    while i > 0 and j > 0: #looks for match/mismatch/gap starting from bottom right of matrix
      if pointerMatrix[i,j] in [2, 5, 6, 9]:
          #appends match/mismatch then moves to the cell diagonally up and to the left
          queryAlign.append(querySequence[j-1])
          subjectAlign.append(subjectSequence[i-1])
          queryAlign.append(qtemp)
          subjectAlign.append(stemp)
          qtemp = ''
          stemp = ''
          i -= 1
          j -= 1
      elif pointerMatrix[i,j] in [4, 6, 7, 9]:
          #appends gap and accompanying nucleotide, then moves to the cell to the left
          queryAlign.append(querySequence[j-1])
          stemp = '-'
          j -= 1
      elif pointerMatrix[i,j] in [3, 5, 7, 9]:
          #appends gap and accompanying nucleotide, then moves to the cell above
          qtemp = '-'
          subjectAlign.append(subjectSequence[i-1])
          i -= 1
      elif pointerMatrix[i,j] == 0:
        queryAlign.append(querySequence[j-1])
        subjectAlign.append(subjectSequence[i-1])
        break
    # Reverses the strings
    subjectAlign = ''.join(subjectAlign)[::-1]
    queryAlign= ''.join(queryAlign)[::-1]
    '\n'.join([subjectAlign, queryAlign])

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

  def __call__(self, subjectSequence: str|list[str], querySequence: str|list[str])->tuple[int,int]: 
      if not numpy:
          raise ImportError('Please pip install numpy!')
      
      if isinstance(subjectSequence, str):
          querySequence,subjectSequence = map(lambda x: x.upper(), [querySequence,subjectSequence])
          querySequence = [x for x in querySequence]
          subjectSequence = [x for x in subjectSequence]

      subjectSequence, querySequence = frontWhiteSpace(subjectSequence, querySequence) 
      #matrix initialisation
      self.alignment_score = numpy.zeros((len(subjectSequence),len(querySequence))) 
      #pointer to trace optimal alignment
      self.pointer = numpy.zeros((len(subjectSequence)+1, len(querySequence)+1)) 
      self.pointer[:,0] = 3
      self.pointer[0,:] = 4

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
  
          if match == tmax: #matrix for traceback based on results from scoring matrix
            self.pointer[i+1,j+1] += 2
          if ugap == tmax:
            self.pointer[i+1,j+1] += 3
          if lgap == tmax:
            self.pointer[i+1,j+1] += 4
  
      return self.alignment_score, self.pointer

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

def frontWhiteSpace(querySequence: list[str], subjectSequence: list[str])->tuple[list[str],list[str]]: 
    #adds leading white space so that matrix works
    subjectSequence = rjustlist(subjectSequence)
    querySequence = rjustlist(querySequence)

    return subjectSequence, querySequence

needlemanWunsch = needleman_wunsch()
watermanSmithBeyer = waterman_smith_beyer()
smithWaterman = smith_waterman()
levenshteinDist = levenshtein()
