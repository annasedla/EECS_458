import numpy as np


# TODO check if the inputs match the signs expected


class Align:

    def __init__(self, alignment, match, mismatch, indel, string_one, string_two):
        self.alignment = alignment
        self.path = 0  # TODO remember to reset path
        self.match = match
        self.mismatch = mismatch
        self.indel = indel
        self.string_one = string_one.upper()
        self.string_two = string_two.upper()

    def align_sequences(self):
        # Store length of two sequences in temp variables
        n = len(self.string_one)
        m = len(self.string_two)

        # Generate matrix of zeros to store scores at each iteration
        score_matrix = np.zeros((len(self.string_two) + 1, len(self.string_one) + 1), np.int)

        if self.alignment == 'g':
            self.needleman_wunsch(n, m, score_matrix)
        elif self.alignment == 'l':
            self.smith_waterman(n, m, score_matrix)
        else:
            print('not a valid alignment')

    # Global alignment algo
    def needleman_wunsch(self, n, m, score_matrix):

        # Calculate score table
        # Fill out first column
        for i in range(0, m + 1):
            score_matrix[i][0] = self.indel * i

        # Fill out first row
        for j in range(0, n + 1):
            score_matrix[0][j] = self.indel * j

        # Fill out all other values in the score matrix
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                # Calculate the max score by checking the top, left, and diagonals
                match = score_matrix[i - 1][j - 1] + self.compute_score(self.string_one[j - 1], self.string_two[i - 1])
                delete = score_matrix[i - 1][j] + self.indel
                insert = score_matrix[i][j - 1] + self.indel

                # Record the maximum score from the three scenarios
                score_matrix[i][j] = max(match, delete, insert)

        print(score_matrix)
        self.num_optimal_solution_nw(len(self.string_two), len(self.string_one), score_matrix)
        print('Num optimal', self.path)
        return self.backtrack(score_matrix)

    # Local alignemnt algo
    def smith_waterman(self, n, m, score_matrix):

        # Fill out all other values in the score matrix
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                # Calculate the max score by checking the top, left, and diagonals
                match = score_matrix[i - 1][j - 1] + self.compute_score(self.string_one[j - 1], self.string_two[i - 1])
                delete = score_matrix[i - 1][j] + self.indel
                insert = score_matrix[i][j - 1] + self.indel

                # Record the maximum score from the three scenarios
                score_matrix[i][j] = max(match, delete, insert)

        print(score_matrix)

        num_indices = np.where(score_matrix == score_matrix.max())

        for i in range(0, len(num_indices[0])):
            self.num_optimal_solution_sw(num_indices[0][i], num_indices[1][i], score_matrix)

            print('Num optimal', self.path)

        print('Final num optimal', self.path)

        return self.backtrack(score_matrix)

    def num_optimal_solution_nw(self, i, j, score_matrix):
        if i == 0 and j == 0:
            self.path = self.path + 1

        elif i == 0:
            self.num_optimal_solution_nw(i, j - 1, score_matrix)

        elif j == 0:
            self.num_optimal_solution_nw(i - 1, j, score_matrix)

        else:
            if (self.string_one[j-1] == self.string_two[i-1]) and (score_matrix[i][j] == score_matrix[i - 1][j - 1] + self.match):
                self.num_optimal_solution_nw(i - 1, j - 1, score_matrix)

            if (self.string_one[j-1] != self.string_two[i-1]) and (score_matrix[i][j] == score_matrix[i - 1][j - 1] + self.mismatch):
                self.num_optimal_solution_nw(i - 1, j - 1, score_matrix)

            if score_matrix[i][j] == score_matrix[i][j-1] + self.indel:
                self.num_optimal_solution_nw(i, j-1, score_matrix)

            if score_matrix[i][j] == score_matrix[i-1][j] + self.indel:
                self.num_optimal_solution_nw(i-1, j, score_matrix)

    def num_optimal_solution_sw(self, i, j, score_matrix):
        if i == 0 or j == 0:
            self.path = self.path + 1

        else:
            if (self.string_one[j - 1] == self.string_two[i - 1]) and (
                    score_matrix[i][j] == score_matrix[i - 1][j - 1] + self.match):
                self.num_optimal_solution_sw(i - 1, j - 1, score_matrix)

            if (self.string_one[j - 1] != self.string_two[i - 1]) and (
                    score_matrix[i][j] == score_matrix[i - 1][j - 1] + self.mismatch):
                self.num_optimal_solution_sw(i - 1, j - 1, score_matrix)

            if score_matrix[i][j] == score_matrix[i][j - 1] + self.indel:
                self.num_optimal_solution_sw(i, j - 1, score_matrix)

            if score_matrix[i][j] == score_matrix[i - 1][j] + self.indel:
                self.num_optimal_solution_sw(i - 1, j, score_matrix)

    # def backtrack_nw(self, score_matrix):
    #
    #     # Traceback and compute the alignment
    #     # Create variables to store alignment
    #     alignment_one = ""
    #     alignment_two = ""
    #
    #     # Start from the bottom right cell in matrix
    #     i = (len(self.string_two))
    #     j = len(self.string_one)
    #
    #     optimal_score = score_matrix[i][j]
    #
    #     while i > 0 and j > 0:
    #         score_current = score_matrix[i][j]
    #         score_diagonal = score_matrix[i - 1][j - 1]
    #         score_up = score_matrix[i][j - 1]
    #         score_left = score_matrix[i - 1][j]
    #
    #         # Check to figure out which cell the current score was calculated from,
    #         # then update i and j to correspond to that cell.
    #
    #         if score_current == score_diagonal + self.compute_score(self.string_one[j - 1], self.string_two[i - 1]):
    #             alignment_one += self.string_one[j - 1]
    #             alignment_two += self.string_two[i - 1]
    #             i -= 1
    #             j -= 1
    #
    #         elif score_current == score_up + self.indel:
    #             alignment_one += self.string_one[j - 1]
    #             alignment_two += '-'
    #             j -= 1
    #
    #         elif score_current == score_left + self.indel:
    #             alignment_one += '-'
    #             alignment_two += self.string_two[i - 1]
    #             i -= 1
    #
    #     # Finish tracing up to the top left cell
    #     while j > 0:
    #         alignment_one += self.string_one[j - 1]
    #         alignment_two += '-'
    #         j -= 1
    #     while i > 0:
    #         alignment_one += '-'
    #         alignment_two += self.string_two[i - 1]
    #         i -= 1
    #
    #     # reverse the sequences
    #     alignment_one = alignment_one[::-1]
    #     alignment_two = alignment_two[::-1]
    #
    #     return alignment_one, alignment_two, optimal_score

    def backtrack(self, score_matrix):

        # Traceback and compute the alignment
        # Create variables to store alignment
        alignment_one = ""
        alignment_two = ""

        if self.alignment == 'g':
            # Start from the bottom right cell in matrix
            i, j = (len(self.string_two)), len(self.string_one)
            optimal_score = score_matrix[i][j]

        elif self.alignment == 'l':
            # Start from the highest element in the matrix
            i, j = np.unravel_index(score_matrix.argmax(), score_matrix.shape)
            optimal_score = score_matrix.max()

        else:
            i, j = 0, 0
            optimal_score = ''
            print('Not a valid alignment symbol')

        while i > 0 and j > 0:
            score_current = score_matrix[i][j]
            score_diagonal = score_matrix[i - 1][j - 1]
            score_up = score_matrix[i][j - 1]
            score_left = score_matrix[i - 1][j]

            # Check to figure out which cell the current score was calculated from,
            # then update i and j to correspond to that cell.

            if score_current == score_diagonal + self.compute_score(self.string_one[j - 1], self.string_two[i - 1]):
                alignment_one += self.string_one[j - 1]
                alignment_two += self.string_two[i - 1]
                i -= 1
                j -= 1

            elif score_current == score_up + self.indel:
                alignment_one += self.string_one[j - 1]
                alignment_two += '-'
                j -= 1

            elif score_current == score_left + self.indel:
                alignment_one += '-'
                alignment_two += self.string_two[i - 1]
                i -= 1

        # reverse the sequences
        alignment_one = alignment_one[::-1]
        alignment_two = alignment_two[::-1]

        if self.alignment == 'g':
            # Finish tracing up to the top left cell
            while j > 0:
                alignment_one += self.string_one[j - 1]
                alignment_two += '-'
                j -= 1
            while i > 0:
                alignment_one += '-'
                alignment_two += self.string_two[i - 1]
                i -= 1

        print(alignment_one)
        print(alignment_two)
        print(optimal_score)

    # A function for determining the score between any two bases in alignment
    def compute_score(self, a, b):
        if a == b:
            return self.match
        elif a == '-' or b == '-':
            return self.indel
        else:
            return self.mismatch


def main():
    """
    Output optimal score
    number of solutions that achieve optimal score
    actual alignment
    """

    print("Results:")

    #align = Align('g', 1, -1, -1, 'ATTACA', 'ATGCT')
    align = Align('l', 1, -1, -1, 'catcat', 'cat')
    # output1, output2, optimal_score = align.align_sequences()
    output1, output2, optimal_score = align.align_sequences()

    print('Optimal Score:', optimal_score)
    print('Num Optimal Solutions:')
    print('Actual Alignments:', "\n" + output1 + "\n" + output2)


if __name__ == "__main__":
    main()
