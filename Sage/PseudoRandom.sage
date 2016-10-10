#Blum Blum Shub RNG
""" 
Initializes a Blum-Blum-Shub RNG State.
A BBS-RNG State is a list with two elements:
[N, X] 
N is a 2*bitlen modulus (product of two primes) 
X is the current state of the PRNG.

INPUT:
    bitlen - the bit length of each of the prime factors of n
    seed - a large random integer to start out the prng
OUTPUT:
    state - a BBS-RNG internal state
"""
def BlumBlumShub_Initialize(bitlen, seed): 
# note that this is not the most cryptographically secure
# way to generate primes, because we do not know how the
# internal sage random_prime function works.
  p = 3; 
  while (p < 2^(bitlen-1)) or (3 != (p % 4)): 
    p = random_prime(2^bitlen);
  q = 3; 
  while (q < 2^(bitlen-1)) or (3 != (q % 4)): 
    q = random_prime(2^bitlen);
  N = p*q;
  X = (seed^2 % N)
  state = [N, X]
  return state;

""" Blum-Blum-Shum random number generation function.
INPUT:
  num_bits - the number of bits (iterations) to 
  generate with this RNG.
  state - an internal state of the BBS-RNG (a list [N, X].)
OUTPUT:
  random_bits - a num_bits length list of random bits.
"""
def BlumBlumShub_Generate(num_bits, state): 
  random_bits = [];
  N = state[0] 
  X = state[1]
  for j in xrange(num_bits):
    X = X^2 % N 
    random_bits.append(X % 2)
  # update the internal state 
  state[1] = X;
  return random_bits;

#Linear Congruential RNG
""" 
This functional initializes a linear congruential RNG state.
This state is a list of four integers: [a, c, m, X]
a,c,m are the parameters of the linear congruential instantiation X is the current state of the PRNG.
INPUT:
  a - The coefficient 
  c - The offset 
  m - The modulus 
  X0 - The initial state
OUTPUT:
  state - The initial internal state of the RNG
"""
def LinearCongruential_Initialize(a, c, m, X0):
  return [a,c,m,X0]

""" 
Generates a single linear congruential RNG output and updates the state.
INPUT:
  state - an internal RNG state.
OUTPUT:
  X - a single output of the linear congruential RNG.
"""
def LinearCongruential_Generate(state):
  a = state[0] 
  c = state[1] 
  m = state[2] 
  X = state[3] 
  X_next = (a*X + c) % m 
  state[3] = X_next 
  return X_next