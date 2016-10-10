
#These structures are the underlying
#Galois Field and corresponding Vector Space
#of the field used in the SAES algorithm
#These structures allow us to easily compute with these fields.

#
F = GF(2);
L.<a> = GF(2^4);
V = L.vector_space();
VF8 = VectorSpace(F, 8);
 
#The MixColumns and its Inverse matrices are stored
#as 2x2 matrices with elements in GF(2^4) (as are state matrices.)
#The MixColumns operation (and its inverse) are performed by
#matrix multiplication.
#
MixColumns_matrix = Matrix(L, [[1,a^2],[a^2,1]]);
InverseMixColumns_matrix = MixColumns_matrix.inverse();
SBox_matrix = Matrix(L,
[
[ 		1 + a^3, 			a^2,			a + a^3,  1 + a + a^3],
[ 1 + a^2 + a^3,	  		  1,				a^3,	  1 + a^2],
[ 		a + a^2,			  0,			      a,		1 + a],
[ 	  a^2 + a^3,  a + a^2 + a^3,  1 + a + a^2 + a^3,  1 + a + a^2]
]);


InverseSBox_matrix = Matrix(L,
[
[	a + a^3,	1 + a^2,	  1 + a^3,		 1 + a + a^3],
[		  1,1 + a + a^2,		  a^3, 1 + a + a^2 + a^3],
[	a + a^2,		  0,			a,			   1 + a],
[ a^2 + a^3,		a^2,1 + a^2 + a^3,	   a + a^2 + a^3]
]);

RCON = [
VF8([F(0), F(0), F(0), F(0), F(0), F(0), F(0), F(1)]),
VF8([F(0), F(0), F(0), F(0), F(1), F(1), F(0), F(0)])
]; 

def SAES_ToStateMatrix(block):
	#Converts a bit list into an SAES state matrix.
	B = block;
	# form the plaintext block into a matrix of GF(2^n)	elements
	S00 = L(V([B[0], B[1], B[2], B[3]]));
	S01 = L(V([B[4], B[5], B[6], B[7]]));
	S10 = L(V([B[8], B[9], B[10], B[11]]));
	S11 = L(V([B[12], B[13], B[14], B[15]]));
	state_matrix = Matrix(L, [[S00,S01],[S10,S11]]);
	return state_matrix;

#Converts an SAES State Matrix to a bit list.
def SAES_FromStateMatrix(State_Matrix):
	output = [];
	# convert State Matrix back into bit list
	for r in xrange(2):
		for c in xrange(2):
			v = V(State_Matrix[r,c]);
			for j in xrange(4):
				output.append(Integer(v[j]));
	return output;

def SAES_AddRoundKey(state_matrix, K):
	"""
	Adds a round key to an SAES state matrix.
	"""
	K_matrix = SAES_ToStateMatrix(K);
	next_state_matrix = K_matrix + state_matrix;
	return next_state_matrix;

#Performs the Mix Columns operation.
def SAES_MixColumns(state_matrix):
	next_state_matrix = MixColumns_matrix*state_matrix;
	return next_state_matrix;

#Performs the Inverse Mix Columns operation.
def SAES_InverseMixColumns(state_matrix):
	next_state_matrix = InverseMixColumns_matrix*state_matrix;
	return next_state_matrix;

#Performs the Shift Row operation.
def SAES_ShiftRow(state_matrix):
	M = state_matrix;
	next_state_matrix = Matrix(L, [
	[M[0,0], M[0,1]],
	[M[1,1], M[1,0]]
	]);
	return next_state_matrix;

"""
Performs the SAES SBox look up in the SBox matrix
(lookup table.)
"""
def SAES_SBox(nibble):
	v = nibble._vector_();
	c = Integer(v[0]) + 2*Integer(v[1]);
	r = Integer(v[2]) + 2*Integer(v[3]);
	return SBox_matrix[r,c];

"""
Performs the SAES SBox on each element of an SAES state
matrix.
"""
def SAES_NibbleSubstitution(state_matrix):
	M = state_matrix;
	next_state_matrix = Matrix(L,
	[ [ SAES_SBox(M[0,0]), SAES_SBox(M[0,1])],
	[ SAES_SBox(M[1,0]), SAES_SBox(M[1,1])] ]);
	return next_state_matrix;

"""
Performs the SAES Inverse SBox look up in the SBox
matrix (lookup table.)
"""
def SAES_InvSBox(nibble):
	v = nibble._vector_();
	c = Integer(v[0]) + 2*Integer(v[1]);
	r = Integer(v[2]) + 2*Integer(v[3]);
	return InverseSBox_matrix[r,c];

"""
Performs the SAES Inverse SBox on each element of an
SAES state matrix.
"""
def SAES_InvNibbleSub(state_matrix):
	M = state_matrix;
	next_state_matrix = Matrix(L,
	[ [ SAES_InvSBox(M[0,0]), SAES_InvSBox(M[0,1])],
	[ SAES_InvSBox(M[1,0]), SAES_InvSBox(M[1,1])] ]);
	return next_state_matrix;

def RotNib(w):
	"""
	Splits an 8 bit list into two elements of GF(2^4)
	"""
	N_0 = L(V([w[j] for j in xrange(4)]));
	N_1 = L(V([w[j] for j in xrange(4,8)]));
	return (N_1, N_0);

def SAES_g(w, i):
	"""
	Performs the SAES g function on the 8 bit list w.
	"""
	(N0, N1) = RotNib(w);
	N0 = V(SAES_SBox(N0));
	N1 = V(SAES_SBox(N1));
	temp1 = VF8( [ N0[0], N0[1], N0[2], N0[3],
	N1[0], N1[1], N1[2], N1[3] ] );
	output = temp1 + RCON[i];
	return output;

def SAES_KeyExpansion(K):
	"""
	Expands an SAES key into two round keys.
	"""
	w0 = VF8([K[j] for j in xrange(8)]);
	w1 = VF8([K[j] for j in xrange(8,16)]);
	w2 = w0 + SAES_g(w1, 0);
	w3 = w1 + w2;
	w4 = w2 + SAES_g(w3, 1);
	w5 = w3 + w4;
	K0 = [w0[j] for j in xrange(8)];
	K0.extend([w1[j] for j in xrange(8)]);
	K1 = [w2[j] for j in xrange(8)];
	K1.extend([w3[j] for j in xrange(8)]);
	K2 = [w4[j] for j in xrange(8)];
	K2.extend([w4[j] for j in xrange(8)]);
	return (K0, K1, K2);

#
# Encrypts one plaintext block with key K
#
def SAES_Encrypt(plaintext, K):
	"""
	Performs a SAES encryption on a single plaintext
	block.
	(Both block and key passed as bit lists.)
	"""
	# get the key schedule
	(K0, K1, K2) = SAES_KeyExpansion(K);
	state_matrix0 = SAES_ToStateMatrix(plaintext);
	state_matrix1 = SAES_AddRoundKey(state_matrix0, K0);
	state_matrix2 = SAES_NibbleSubstitution(state_matrix1);
	state_matrix3 = SAES_ShiftRow(state_matrix2);
	state_matrix4 = SAES_MixColumns(state_matrix3);
	state_matrix5 = SAES_AddRoundKey(state_matrix4, K1);
	state_matrix6 = SAES_NibbleSubstitution(state_matrix5);
	state_matrix7 = SAES_ShiftRow(state_matrix6);
	state_matrix8 = SAES_AddRoundKey(state_matrix7, K2);
	output = SAES_FromStateMatrix(state_matrix8);
	return output;

#
# Decrypts one ciphertext block with key K
#
def SAES_Decrypt(ciphertext, K):
	"""
	Performs a single SAES decryption operation on a
	ciphertext block.
	(Both block and key passed as bit lists.)
	"""
	# perform key expansion
	(K0, K1, K2) = SAES_KeyExpansion(K);
	# form the ciphertext block into a matrix of GF(2^n) elements
	state_matrix0 = SAES_ToStateMatrix(ciphertext);
	state_matrix1 = SAES_AddRoundKey(state_matrix0, K2);
	state_matrix2 = SAES_ShiftRow(state_matrix1);
	state_matrix3 = SAES_InvNibbleSub(state_matrix2);
	state_matrix4 = SAES_AddRoundKey(state_matrix3, K1);
	state_matrix5 = SAES_InverseMixColumns(state_matrix4);
	state_matrix6 = SAES_ShiftRow(state_matrix5);
	state_matrix7 = SAES_InvNibbleSub(state_matrix6);
	state_matrix8 = SAES_AddRoundKey(state_matrix7, K0);
	output = SAES_FromStateMatrix(state_matrix8);
	return output;