
alph = { 'A', 'C', 'G', 'U', 'N'}

class RNA_Fold():
  """
  A class to solve the optimal secondary structure for a given RNA sequence
  """

  def __init__(self, s=""):
    self.rna = s
    self.n = len(s)
    self.F = [[ '!' for i in range(self.n) ] for j in range(self.n) ]


  def fold(self,s):
    """
    Given a string s representing an RNA sequence
    determines the optimal secondary structure via dynamic programming

    fills a dp matrix with optimal subproblem values

    returns the optimal score, and once backtracking is implemented will return the optimal fold
    """
    self.rna = s
    self.n = len(s)
    self.F = [[ '!' for i in range(self.n) ] for j in range(self.n) ]

    return self.f_iter()


  def f(self,i,j):
    """
    main recursive dp call, f(i,j) represents the value of the optimal secondary structure 
    of the rna sequence between positions i and j 
    """
    #print("i: " + str(i) + " j: " + str(j))

    #Memoize results
    if(self.F[i][j] != '!'):
      return self.F[i][j]

    #Base Case 1) no bases can be aligned
    if( i >= j or i == j-1):
      self.F[i][j] = 0
      return self.F[i][j]


    case1 = self.f(i+1,j)
    case2 = self.f(i,j-1)
    case3 = 0
    case4 = 0

    if(self.complements( self.rna[i], self.rna[j] )):
      case3 = self.f(i+1,j-1) + 1

    for k in range(j-i):
      s = self.f(i,i+k) + self.f(i+k+1,j)
      if(s > case4):
        case4 = s

    self.F[i][j] = max(case1, case2, case3, case4)
    return self.F[i][j]

  def f_iter(self):
    #base
    for i in range(self.n):
      self.F[i][i] = 0
      self.F[i][i-1] = 0

    for g in range(1,self.n):
      for i in range(self.n-g):
        j = i + g
        case1 = self.F[i+1][j]
        case2 = self.F[i][j-1]
        case3 = 0
        case4 = 0
        if self.complements(self.rna[i], self.rna[j] ):
          case3 = self.F[i+1][j-1] + 1 
        for k in range(j-i):
          s = self.F[i][i+k] + self.F[i+k+1][j]
          if( s > case4):
            case4 = s

        self.F[i][j] = max(case1, case2, case3, case4)

    return self.F[0][self.n-1]

  def complements(self,x, y):
    """
    Given characters x and y in alph, returns true is x is the complement of y
    """
    #TODO should this thrown an error? probably
    if(x not in alph or y not in alph): return False
    if (x == 'N' or y == 'N'): return True
    if(x == 'A' and y == 'U'): return True
    if(x == 'U' and y == 'A'): return True
    if(x == 'G' and y == 'C'): return True
    if(x == 'C' and y == 'G'): return True
    return False

  def print_F(self):

    space =  len(str(max(self.F[0])))

    print(space)

    line = " "*space
    for c in self.rna:
      line += c + " "*space
    print(line)
    line = ""
    for j in range(self.n):
      line += self.rna[j] + " "*space
      for i in range(self.n):
        line += str(self.F[j][i]) + " "*space
      print(line)
      line = ""
    

