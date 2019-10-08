import numpy as np


class Align:

    def __init__(self):
        """
        Global variables
        """
        self.alignment = ''
        self.path = 0
        self.match = 1
        self.mismatch = -1
        self.indel = -1
        self.string_one = ''
        self.string_two = ''

    def align_sequences(self, alignment, string_one, string_two, path=0, match=1, mismatch=-1, indel=-1):
        """
        Method initially called
        :param alignment: l or g
        :param string_one: x
        :param string_two: y
        :param path: total bumber of optimal paths
        :param match:
        :param mismatch:
        :param indel:
        :return: nothing, but prints to the console
        """

        # Set all the global variables
        self.alignment = alignment
        self.string_one = string_one.upper()
        self.string_two = string_two.upper()
        self.path = path
        self.match = match
        self.mismatch = mismatch
        self.indel = indel

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
        """
        Global alignment
        :param n: rows
        :param m: columns
        :param score_matrix: initial score matrix
        :return: prints all the results
        """

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

        """
        OUTPUT
        """
        print("RESULTS:")
        print()
        print('~~Global alignment~~')
        print("Score Matrix:")
        print(score_matrix)

        self.num_optimal_solution_nw(len(self.string_two), len(self.string_one), score_matrix)

        print('Number of optimal solutions:', self.path)
        self.backtrack_nw(score_matrix)

    # Local alignemnt algo
    def smith_waterman(self, n, m, score_matrix):
        """
        Local alignment
        :param n: rows
        :param m: columns
        :param score_matrix: initial score matrix
        :return: prints all the results
        """

        # Fill out all other values in the score matrix
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                # Calculate the max score by checking the top, left, and diagonals
                match = score_matrix[i - 1][j - 1] + self.compute_score(self.string_one[j - 1], self.string_two[i - 1])
                delete = score_matrix[i - 1][j] + self.indel
                insert = score_matrix[i][j - 1] + self.indel

                # Record the maximum score from the three scenarios
                score_matrix[i][j] = max(match, delete, insert, 0)

        num_indices = np.where(score_matrix == score_matrix.max())

        for i in range(0, len(num_indices[0])):
            self.num_optimal_solution_sw(num_indices[0][i], num_indices[1][i], score_matrix)

        """
        OUTPUT
        """
        print("RESULTS:")
        print()
        print('~~Local alignment~~')
        print("Score Matrix:")
        print(score_matrix)

        print('Number of optimal solutions:', self.path)
        self.backtrack_sw(score_matrix)

    def num_optimal_solution_nw(self, i, j, score_matrix):
        """
        Computes the number of total optimal paths for global alignment
        :param i: initially length of string two
        :param j: initially length of string one
        :param score_matrix: finalized score matrix
        :return: number of optimal paths
        """
        if i == 0 and j == 0:
            self.path = self.path + 1

        elif i == 0:
            self.num_optimal_solution_nw(i, j - 1, score_matrix)

        elif j == 0:
            self.num_optimal_solution_nw(i - 1, j, score_matrix)

        else:
            if (self.string_one[j-1] == self.string_two[i-1]) and (score_matrix[i][j] == score_matrix[i - 1][j - 1] +
                                                                   self.match):
                self.num_optimal_solution_nw(i - 1, j - 1, score_matrix)

            if (self.string_one[j-1] != self.string_two[i-1]) and (score_matrix[i][j] == score_matrix[i - 1][j - 1] +
                                                                   self.mismatch):
                self.num_optimal_solution_nw(i - 1, j - 1, score_matrix)

            if score_matrix[i][j] == score_matrix[i][j-1] + self.indel:
                self.num_optimal_solution_nw(i, j-1, score_matrix)

            if score_matrix[i][j] == score_matrix[i-1][j] + self.indel:
                self.num_optimal_solution_nw(i-1, j, score_matrix)

    def num_optimal_solution_sw(self, i, j, score_matrix):
        """
        Computes the number of total optimal paths for smith waterman
        :param i: initially length of string two
        :param j: initially length of string one
        :param score_matrix: finalized score matrix
        :return: number of optimal paths
        """

        if score_matrix[i][j] == 0:
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

    def backtrack_nw(self, score_matrix):
        """
        Prints out the final string alignments
        :param score_matrix: Matrix filled with values
        :return: finalized string alignments
        """

        # Traceback and compute the alignment
        # Create variables to store alignment
        alignment_one = ""
        alignment_two = ""

        # Start from the bottom right cell in matrix
        i, j = (len(self.string_two)), len(self.string_one)
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

            elif score_current == score_up + self.indel:
                alignment_one += self.string_one[j - 1]
                alignment_two += '-'
                j -= 1

            elif score_current == score_left + self.indel:
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

        print('String one:', alignment_one)
        print('String two:', alignment_two)
        print('Optimal score:', optimal_score)

    def backtrack_sw(self, score_matrix):
        """
        Prints out the final string alignments
        :param score_matrix: Matrix filled with values
        :return: finalized string alignments
        """

        # Traceback and compute the alignment
        # Create variables to store alignment
        alignment_one = ""
        alignment_two = ""

        # Start from the highest element in the matrix
        i, j = np.unravel_index(score_matrix.argmax(), score_matrix.shape)
        optimal_score = score_matrix.max()

        while score_matrix[i][j] != 0:
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

        print('String one:', alignment_one)
        print('String two:', alignment_two)
        print('Optimal score:', optimal_score)

    def compute_score(self, a, b):
        """
        A function for determining the score between any two bases in alignment
        :param a: character one
        :param b: character two
        :return: whether or not they are the same
        """
        if a == b:
            return self.match
        elif a == '-' or b == '-':
            return self.indel
        else:
            return self.mismatch


def main():

    align = Align()  # make an instance of a class

    alignment = ''
    match = 0
    mismatch = 0
    indel = 0
    string_one = ''
    string_two = ''

    filepath = 'test_input.txt'  # place the correct filepath here

    with open(filepath) as fp:
        line = fp.readline()
        cnt = 1
        while line:
            if cnt == 1:
                alignment = str(line.strip())

            if cnt == 2:
                current_line = line.strip().split(',')
                match = int(current_line[0])
                mismatch = int(current_line[1])
                indel = int(current_line[2])

            if cnt == 3:
                string_one = str(line.strip())

            if cnt == 4:
                string_two = str(line.strip())

            line = fp.readline()
            cnt += 1

    align.align_sequences(alignment, string_one, string_two, match=match, mismatch=mismatch, indel=indel)


if __name__ == "__main__":
    main()
