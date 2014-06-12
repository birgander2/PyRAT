def primefactors(n):
    """
    Calculates a factorisation into prime numbers
    """
    result = []
    for i in range(2,n+1):
        s = 0
        while n%i == 0: 
            n /= i
            s += 1
        result.extend([i]*s) 
        if n==1:
            return result
