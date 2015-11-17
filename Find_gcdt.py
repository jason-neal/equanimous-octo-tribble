#!/usr/bin/env python
#Find_gcdt.py
# Find a common step size dt of a timeseries.
import numpy as np
import decimal

def gcdt(ts, dp=False):
    """Find the greatest common step size dt.
    inputs:
        ts - the timeseries for which to find the step dt.
        dp - the number of decimal places accuracy required
            if not provided it will calcualte the decimal precision in 
            the values
    output: 
        dt - timestep dt that is common to all values.     
    """   
    if not dp:
        numdp = np.zeros_like(ts)
        for i, val in enumerate(ts):
            numdp[i] = abs(decimal.Decimal(str(val)).as_tuple().exponent)   
        dp = np.max(numdp)
        print("The max number of dp found was ", dp)
    
    multiplier = float(10 ** dp)
    #print(ts)
    new_list = [int(val * multiplier) for val in ts]
    #print(new_list)
    count = 0
    while len(new_list) > 1:
    	print("Counter = ", count)
        old_list = new_list
        new_list = []
    # Iterate common divisors until only one left or gcd is 1.
        for i in range(len(old_list)-1):
            cd = gcd(old_list[i], old_list[i+1])
            if cd == 1: #gcd is 1
                #print("gcd found was 1")
                return cd / multiplier
            else:
                new_list.append(cd)
        new_list = list(set(new_list))  # unique values only
        print("unique gcd", new_list)
                
        count += 1
    dt = new_list[0]  # last common divisor   
    print(multiplier)
    return dt/multiplier

def gcd(a, b):
    """Calculate the Greatest Common Divisor of a and b.

    Unless b==0, the result will have the same sign as b (so that when
    b is divided by it, the result comes out positive).
    """
    while b:
        a, b = b, a%b
    return a


if __name__ == "__main__":
 
    ts = [0.0122, 0.1072, 0.30052, 0.7272, 1.28]
    # ts = [random.randint(1, 2035)/3. for i in range(10)]
    
    dt1 = gcdt(ts)
    dt2 = gcdt(ts, 3)
    dt3 = gcdt(ts, 7)
    #print("the answer should be 0.05")
    print("Common dt found floating =", dt1)
    print("Common dt found at 3dp =", dt2)
    print("Common dt found at 7dp =", dt3)