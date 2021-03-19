
def generator(func):
    def inner_func(*args):
        number1 = args[0]
        number2 = args[1]
        return func(number1, number2)
    return inner_func

@generator
def add(*args):
    print(args[0])
    print(args[1])
    return "Finished"

def main():
    print("running program")
    print(add(1,2,1))

if __name__ == '__main__':
    main()