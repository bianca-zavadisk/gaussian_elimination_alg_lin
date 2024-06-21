class DimensionError(Exception):
    def __init__(self, message, invalid_value):
        super().__init__(message)
        self.invalid_value = invalid_value

class Vector():
    def __init__(self, dim: int, coord: list = []):
        self.dim = dim
        self.coord = coord
        
        if self.coord == []:
            for k in range(dim):
                self.coord.append([0]*dim)
        
    def __mul__(self, scalar):
        scaled_vector = Vector(self.dim, self.coord)
        for k in range(self.dim):
            scaled_vector.coord[k] *= scalar
        
        return scaled_vector
    
    def insert(self, entry: int, value: float):
        new_vector = Vector(self.dim, self.coord)
        new_vector.coord[entry-1] = value
        
        return new_vector
    
    def __str__(self) -> str:
        return self.coord
    
class Matrix():
    def __init__(self, i: int, j: int, entries: list = []) -> None:
        self.rows = i
        self.columns = j
        self.entries = entries
        
        if self.entries == []:
            self.entries = [[0.0] * self.columns for _ in range(self.rows)]
            
    def insert(self, line: int, column: int, value: float):
        new_matrix = Matrix(self.rows, self.columns, self.entries)
        new_matrix.entries[line-1][column-1] = value
        
        return new_matrix
        
    def __add__(self, other_matrix):
        if self.rows != other_matrix.rows or self.columns != other_matrix.columns:
            raise DimensionError('Both matrixes should have the same dimension', other_matrix)
        
        sum_matrix = Matrix(self.rows, self.columns)
        
        for i in range(self.rows):
            for j in range(self.columns):
                sum_matrix.entries[i][j] = self.entries[i][j] + other_matrix.entries[i][j]
        
        return sum_matrix
    
    def __sub__(self, other_matrix):
        if self.rows != other_matrix.rows or self.columns != other_matrix.columns:
            raise DimensionError('Both matrixes should have the same dimension', other_matrix)
        
        sub_matrix = Matrix(self.rows, self.columns)
        
        for i in range(self.rows):
            for j in range(self.columns):
                sub_matrix.entries[i][j] = self.entries[i][j] - other_matrix.entries[i][j]
        
        return sub_matrix
    
    def __mul__(self, element):
        #Multiplication by scalar
        if isinstance(element, int) or isinstance(element, float):
            mul_matrix = Matrix(self.rows, self.columns)
            for i in range(self.rows):
                for j in range(self.columns):
                    mul_matrix.entries[i][j] = self.entries[i][j]*element
        
        #Matrix multiplication
        if isinstance(element, Matrix):
            if self.columns != element.rows:
                raise DimensionError('Matrix column number should match the seccond matrix rows number', element.rows)
            mul_matrix = Matrix(self.rows, element.columns)
            for i in range(self.rows):
                for j in range(element.columns):
                    sum_of_products = 0
                    for k in range(self.rows):
                        # mul_matrix.entries[i][j] += self.entries[i][k]*element.entries[k][j]
                        sum_of_products += self.entries[i][k]*element.entries[k][j]
                    mul_matrix.entries[i][j] = sum_of_products
        
        else:
            raise TypeError('Multiplication not supported', element)
        
        return mul_matrix
    
    def __eq__(self, other_matrix):
        equal = False
        if self.rows == other_matrix.rows and self.columns == other_matrix.columns:
            equal = True
        for i in range(self.rows):
            for j in range(self.columns):
                if self.entries[i][j] != other_matrix.entries[i][j]:
                    equal = False
                    return equal
        
        return equal
    
    def __pow__(self, power: int):
        pow_matrix = self
        
        if power < -1:
            raise ValueError('Power with negative exponents are not definide')
        
        if power == -1:
            raise ValueError('Inverse is not defined')
        
        if power == 0:
            return pow_matrix.identity()
        
        for i in range(power - 1):
            pow_matrix *= self
            
        return pow_matrix
    
    def transpose(self):
        transpose_matrix = Matrix(self.columns, self.rows)
        for i in range(self.rows):
            for j in range(self.columns):
                transpose_matrix.entries[j][i] = self.entries[i][j]
        
        return transpose_matrix
    
    def identity(self):
        identity_matrix = Matrix(self.rows, self.columns)
        if self.rows != self.columns:
            raise DimensionError('Matrix should be square')
        
        for i in range(self.rows):
            identity_matrix.entries[i][i] = 1.0  # Set diagonal elements to 1
            
        return identity_matrix
    
    def __str__(self) -> str:
        str_matrix = ''
        for k in range(self.rows):
            str_matrix += f'{self.entries[k]}\n'
        return str_matrix      

def permutation(dim: int, line_i: int, line_j: int) -> Matrix:
    p_matrix = Matrix(dim, dim, [])    
    id_matrix = Matrix(dim, dim).identity()
    
    for k in range(dim):
        p_matrix.entries[k] = id_matrix.entries[k]
        if k == line_i:
            p_matrix.entries[k] = id_matrix.entries[line_j]
        if k == line_j:
            p_matrix.entries[k] = id_matrix.entries[line_i]

    return p_matrix

def line_scale(dim: int, line_i: int, scale: float, inverse=False) -> Matrix:
    if scale == 0:
        raise ValueError('Operation not possible')
    
    line_scale_matrix = Matrix(dim, dim).identity()
    
    line_scale_matrix.entries[line_i][line_i] *= scale
    
    if inverse:
        line_scale_matrix.entries[line_i][line_i] = 1/(scale)
    
    return line_scale_matrix

def lines_operation(dim: int, line_i: int, line_j: int, scale=-1.0, inverse=False) -> Matrix:
    lines_operation_matrix = Matrix(dim, dim).identity()
    
    lines_operation_matrix.entries[line_i][line_j] = scale
    
    if inverse:
        lines_operation_matrix.entries[line_i][line_j] = (-1)*scale
    
    return lines_operation_matrix

def permut_l(l_matrix, line_i, line_j, till_j):
    for j in range(till_j):
        l_matrix.entries[line_j][j], l_matrix.entries[line_i][j] = l_matrix.entries[line_i][j], l_matrix.entries[line_j][j]
        
    return l_matrix

def lower_triangular(dim: int, previous_l: Matrix, op_matrix: Matrix) -> Matrix:
    l_matrix = op_matrix*previous_l
    
    return l_matrix

def gaussian_elimination(a: list):
    '''
    Performs Gaussian elimination on a n x m matrix A and returns the permutation matrix P, lower triangular matrix L, and upper triangular matrix U.
    
    Inputs:
    a: list - n x m matrix A
    
    Outputs:
    p: Matrix - permutation matrix P
    l: Matrix - lower triangular matrix L
    u: Matrix - upper triangular matriix U
    '''
    
    a = Matrix(len(a), len(a[0]), a) #Define Matrix A as an element of the class Matrix(), in order to use pre-defined operations
    
    n = a.rows
    m = a.columns
    
    p_matrix = Matrix(n, n).identity()
    u_matrix = Matrix(n, m, a.entries)
    l_matrix = lower_triangular(n, Matrix(n,n).identity(), Matrix(n,n).identity())
    
    for i in range(n):
        for j in range(m):
            if i == j:
                if u_matrix.entries[i][j] == 0:
                    pivot_row = j
                    for k in range (i+1, n):
                        if abs(u_matrix.entries[k][j]) > abs(u_matrix.entries[pivot_row][j]):
                            pivot_row = k
                
                    #Row_exchange
                    p_op_matrix = permutation(n, i, pivot_row)
                    p_matrix = p_op_matrix*p_matrix
                    u_matrix = p_op_matrix*u_matrix
                    
                    l_matrix = permut_l(l_matrix, i, pivot_row, i)
        
                pivot = u_matrix.entries[i][j]      
        
        if pivot == 0:
            continue
        
        min_piv = min(n, m) #When A is not square, the number of pivots is the minimun between n, m    

        #Row operations
        for k in range(i+1, min_piv):     
            scale = (u_matrix.entries[k][i] / pivot)
            if scale:
                op_matrix = lines_operation(n, k, i, -scale)
                u_matrix = op_matrix*u_matrix
                l_matrix.entries[k][i] = l_matrix.entries[i][i]*scale
                
    #Row scale to get pivot = 1; After all operations are done
    for i in range(n):
        for j in range(m):
            if i == j:
                pivot = u_matrix.entries[i][j]
                    
        if not pivot: #When pivot = 0, uses the next non-zero element to scale the row
            for k in range(i+1, m):
                if u_matrix.entries[i][k]:
                    pivot = u_matrix.entries[i][k]
                    break
        
        if pivot != 1 and pivot:
            op_matrix = line_scale(n, i, (1/pivot))
            u_matrix = op_matrix*u_matrix
            for l in range(n):
                l_matrix.entries[l][i] *= pivot
                if l_matrix.entries[l][i] == -0.0:
                    l_matrix.entries[l][i] = 0.0
                
            
        for k in range(1, n-i): #Divides linearly dependent rows, when equal
            if u_matrix.entries[i] == u_matrix.entries[i+k]:
                op_matrix = lines_operation(n, i, i+k, -1)
                u_matrix = op_matrix*u_matrix
                l_matrix.entries[i][i+k] = 1
    
    decomp_lu = [p_matrix, l_matrix, u_matrix]
    
    #Print P, L and U in matrix form
    print(f'P:\n{p_matrix} \nL:\n{l_matrix}\nU:\n{u_matrix}')
        
    #Returns a list containing P, L and U in type() = Matrix, maintain operations and properties from the class
    return decomp_lu

    #If you want just a list of lists inputs, it can be replaced by:
    decomp_lu_entries = [p_matrix.entries, l_matrix.entries, u_matrix.entries]
    return decomp_lu_entries