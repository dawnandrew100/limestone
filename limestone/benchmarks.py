from random import randint, choice
from string import ascii_lowercase
from timeit import timeit
import limestone as ls

test_strings1 = ["".join(choice(ascii_lowercase) for _ in range(randint(50,100))) for _ in range(1000)]
test_strings2 = ["".join(choice(ascii_lowercase) for _ in range(randint(50,100))) for _ in range(1000)]

def test_NW():
    nw_list = [ls.needleman_wunsch.distance(a, b) for a, b in zip(test_strings1, test_strings2)]

def main():
    print(f"Needleman Wunsch test: Time ={timeit('test_NW()', globals = globals(), number=1)}")

if __name__ == "__main__":
    main()
