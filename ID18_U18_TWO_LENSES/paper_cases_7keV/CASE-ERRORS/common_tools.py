def shuffle(x, interval=0.1):
    import numpy.random
    xnew = x + (numpy.random.random() - 0.5) * (x * interval)
    if isinstance(x, int):
        return int(xnew)
    else:
        return xnew

print(shuffle(800))
print(shuffle(12.0))