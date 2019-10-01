import numpy as np


class Alignment:

    def __init__(self, string_one, string_two):
        self.string_one = string_one
        self.string_two = string_two

        self.matrix = np.zeros((len(string_one), len(string_two)))

    def populate_matrix(self):
        for i in range(0, len(self.string_one)):
            print('hi')


def main():
    print("hi")

if __name__ == "__main__":
    main()
