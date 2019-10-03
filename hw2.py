import numpy as np

# constants to calculate scores
gap_penalty = -1  # indel
match_award = 1  # match
mismatch_penalty = -1  # mismatch # TODO rewrite this

# TODO check if the inputs match the signs expected


class Align:

    def __init__(self, alignment, match, mismatch, indel, string_one, string_two):
        self.alignment = 'l'
        self.path = 0
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

        if self.alignment == 'l':
            return self.needleman_wunsch(n, m, score_matrix)
        else:
            return 0

    # # Local alignemnt algo
    # def smith_waterman(self, n, m, score_matrix):
    #     # Fill out all other values in the score matrix
    #     for i in range(1, m + 1):
    #         for j in range(1, n + 1):
    #             # Calculate the max score by checking the top, left, and diagonals
    #             match = score_matrix[i - 1][j - 1] + self.compute_score(self.string_one[j - 1], self.string_two[i - 1])
    #             delete = score_matrix[i - 1][j] + gap_penalty
    #             insert = score_matrix[i][j - 1] + gap_penalty
    #
    #             # Record the maximum score from the three scenarios
    #             score_matrix[i][j] = max(match, delete, insert)
    #
    #     print(score_matrix)
    #
    #     b, optimal_score = self.backtrack_sw(score_matrix)
    #
    # def backtrack_sw(self, score_matrix):


    # Global alignment algo
    def needleman_wunsch(self, n, m, score_matrix):

        # # Store length of two sequences in temp variables
        # n = len(self.string_one)
        # m = len(self.string_two)
        #
        # # Generate matrix of zeros to store scores at each iteration
        # score_matrix = np.zeros((len(self.string_two) + 1, len(self.string_one) + 1), np.int)

        # Calculate score table
        # Fill out first column
        for i in range(0, m + 1):
            score_matrix[i][0] = gap_penalty * i

        # Fill out first row
        for j in range(0, n + 1):
            score_matrix[0][j] = gap_penalty * j

        # Fill out all other values in the score matrix
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                # Calculate the max score by checking the top, left, and diagonals
                match = score_matrix[i - 1][j - 1] + self.compute_score(self.string_one[j - 1], self.string_two[i - 1])
                delete = score_matrix[i - 1][j] + gap_penalty
                insert = score_matrix[i][j - 1] + gap_penalty

                # Record the maximum score from the three scenarios
                score_matrix[i][j] = max(match, delete, insert)

        print(score_matrix)
        return self.backtrack_nw(score_matrix)

    def num_optimal_solution(self, i, j, score_matrix, path):
        if i == 1 and j == 1:
            return path + 1

        elif i == 1:
            return self.num_optimal_solution(i, j - 1, score_matrix, path)

        elif j == 1:
            return self.num_optimal_solution(i - 1, j, score_matrix, path)

        else:
            if (self.string_one(i) == self.string_two(i)) and (score_matrix[i][j] == score_matrix[i - 1][j - 1] + match_award):
                return self.num_optimal_solution(i - 1, j - 1, score_matrix, path)

            if (self.string_one(i) != self.string_two(i)) and (score_matrix[i][j] == score_matrix[i - 1][j - 1] + mismatch_penalty):
                return self.num_optimal_solution(i - 1, j - 1, score_matrix, path)

            if score_matrix[i][j] == score_matrix[i][j-1] + gap_penalty:
                return self.num_optimal_solution(i, j-1, score_matrix, path)

            if score_matrix[i][j] == score_matrix[i-1][j] + gap_penalty:
                return self.num_optimal_solution(i-1, j, score_matrix, path)

            return 0

    def backtrack_nw(self, score_matrix):

        # Traceback and compute the alignment
        # Create variables to store alignment
        alignment_one = ""
        alignment_two = ""

        # Start from the bottom right cell in matrix
        i = (len(self.string_two))
        j = len(self.string_one)

        optimal_score = score_matrix[i][j]

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

            elif score_current == score_up + gap_penalty:
                alignment_one += self.string_one[j - 1]
                alignment_two += '-'
                j -= 1

            elif score_current == score_left + gap_penalty:
                alignment_one += '-'
                alignment_two += self.string_two[i - 1]
                i -= 1

        # Finish tracing up to the top left cell
        while j > 0:
            alignment_one += self.string_one[j - 1]
            alignment_two += '-'
            j -= 1
        while i > 0:
            alignment_one += '-'
            alignment_two += self.string_two[i - 1]
            i -= 1

        # reverse the sequences
        alignment_one = alignment_one[::-1]
        alignment_two = alignment_two[::-1]

        return alignment_one, alignment_two, optimal_score

    # A function for determining the score between any two bases in alignment
    def compute_score(self, a, b):
        if a == b:
            return match_award
        elif a == '-' or b == '-':
            print('here')
            return gap_penalty
        else:
            return mismatch_penalty


def main():
    """
    Output optimal score
    number of solutions that achieve optimal score
    actual alignment
    """

    print("Results:")

    #align = Align('g', -1, -1, -1, 'ATTACA', 'ATGCT')
    align = Align('g', -1, -1, -1, 'catcat', 'cat')
    output1, output2, optimal_score = align.align_sequences()

    print('Optimal Score:', optimal_score)
    print('Num Optimal Solutions:')
    print('Actual Alignments:',"\n" + output1 + "\n" + output2)


if __name__ == "__main__":
    main()
