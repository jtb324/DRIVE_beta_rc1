def add_decorator(name):
    def add(func):
        def inner_func(*args, **kwargs):
            print(f"Hello {name}")
            print(f"the analysis is {kwargs['analysis']}")
            return func(*args, **kwargs)
        return inner_func
    return add

@add_decorator("JT")
def add(num1, num2, analysis):
    print(analysis)
    print(num1+num2)

add(3,5,analysis="add")