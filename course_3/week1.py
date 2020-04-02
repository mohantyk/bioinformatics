from math import inf

def dp_change(money, coins):
    '''
    Dynamic programming implementation for calculating minimum coins
    to make change for money.
    '''
    results = [inf]*(money+1)
    results[0] = 0
    for idx in range(1, money+1):
        for coin in coins:
            if coin <= idx:
                result = 1+results[idx-coin]
                if result < results[idx]:
                    results[idx] = result
    return results[money]