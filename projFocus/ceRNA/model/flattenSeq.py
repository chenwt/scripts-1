from collections import Iterable
def flatten(items, ignore_types=(str, bytes)): 
    for x in items:
        if isinstance(x, Iterable) and not isinstance(x, ignore_types): 
            for element in flatten(x):
                yield element
                # yield from flatten(x)
        else:
            yield x
def __test__():
    # items = [1, 2, [3, 4, [5, 6], 7], 8]
    items = {'x':1,'y':[2,3,4]}
    for x in flatten(items): 
        print(x)
def __main__():
    __test__()
