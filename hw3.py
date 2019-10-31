

class Viterbi:

    def __init__(self, observations, states, initial_probabilities, transition_probabilities, emission_probabilities):
        """
        Initialization method
        :param observations: List of observations, aka heads or tails
        :param states: In this case fair or biased
        :param initial_probabilities: probability of picking fair vs biased coin
        :param transition_probabilities: probability that the dealer switches
        :param emission_probabilities: probability of heads or tails for each coin
        """
        self.observations = observations
        self.states = states
        self.initial_probabilities = initial_probabilities
        self.transition_probabilities = transition_probabilities
        self.emission_probabilities = emission_probabilities

        # matrrices
        self.viterbi_matrix = {}
        self.backward_matrix = {}
        self.forward_matrix = {}
        self.viterbi_path = [None] * len(self.observations)
        self.y_probability = 0

    def viterbi(self):
        """
        Creates a Viterbi matrix and fills it up
        :return: N/A
        """
        for state in self.initial_probabilities:
            self.viterbi_matrix[state] = [None] * (len(self.observations))
            self.viterbi_matrix[state][0] = (self.initial_probabilities[state] *
                                             self.emission_probabilities[state][self.observations[0]], "S")

        for i in range(len(self.observations)-1):
            for transition in self.transition_probabilities:
                max = 0
                candidate_state = None
                for state in self.states:
                    previous_probability = self.viterbi_matrix[state][i][0]
                    transition_probability = self.transition_probabilities[state][transition]
                    emission_probability = self.emission_probabilities[transition][self.observations[i + 1]]

                    candidate = previous_probability * transition_probability * emission_probability

                    if candidate > max:
                        max = candidate
                        candidate_state = state

                self.viterbi_matrix[transition][i + 1] = (max, candidate_state)

    def print_viterbi_matrix(self):
        """
        Method to print matrices in a nice way
        :return: N/A
        """
        matrix = [[0 for x in range(len(self.observations))] for y in range(len(self.states))]

        i = 0

        for state in self.states:
            for observation in range(len(self.observations)):
                matrix[i][observation] = round(self.viterbi_matrix[state][observation][0], 6)
            i = i + 1

        print(matrix)

    def compute_viterbi_path(self):
        """
        Computes the viterbi path for a coin used
        :return: N/A
        """
        previous_state = None
        probability = 0

        for state in self.states:
            if self.viterbi_matrix[state][-1][0] > probability:
                probability = self.viterbi_matrix[state][-1][0]
                previous_state = self.viterbi_matrix[state][-1][1]

        for i in range(len(self.observations)):
            self.viterbi_path[-(i+1)] = previous_state
            previous_state = self.viterbi_matrix[previous_state][-(i+1)][1]

        print(self.viterbi_path)

    def calculate_forward_probabilities(self):
        """
        Calculates forward probabilities
        :return: N?A
        """

        y_probability = 0

        for state in self.initial_probabilities:
            self.forward_matrix[state] = [0] * (len(self.observations))
            self.forward_matrix[state][0] = \
                self.initial_probabilities[state]*self.emission_probabilities[state][self.observations[0]]

        for i in range(len(self.observations)-1):
            for state in self.states:
                for transition in self.transition_probabilities:
                    observation = self.observations[i+1]
                    emission_probability = self.emission_probabilities[state][observation]
                    transition_probabilty = self.transition_probabilities[transition][state]
                    previous_probability = self.forward_matrix[transition][i]
                    self.forward_matrix[state][i+1] += \
                        emission_probability * transition_probabilty * previous_probability

        for state in self.states:
            y_probability += self.forward_matrix[state][-1]

        self.y_probability = y_probability

        print(self.forward_matrix)

    def calculate_backward_probabilities(self):
        """
        Calculates backward probabilities
        :return: N/A
        """
        for state in self.states:
            self.backward_matrix[state] = [0] * (len(self.observations))
            self.backward_matrix[state][-1] = 1

        for i in range(len(self.observations)):
            for state in self.states:
                for transition in self.transition_probabilities:
                    observation = self.observations[-(i+1)]
                    emission_probability = self.emission_probabilities[state][observation]
                    transition_probability = self.transition_probabilities[state][transition]
                    previous_probability = self.backward_matrix[transition][-i]
                    self.backward_matrix[transition][-(i+1)] +=\
                        emission_probability * transition_probability * previous_probability

        print(self.backward_matrix)

    def calculate_posterior_probability(self, position, state):
        """
        Computes the posterior based on the forward and backward probabilities
        :param position: what position in the tosses
        :param state: fair vs biased
        :return: N/A
        """
        forward_state = self.forward_matrix[state][position]
        backward_state = self.backward_matrix[state][position]
        posterior = (forward_state*backward_state)/self.y_probability
        print(posterior)


def main():

    '''
    PROBLEM SETUP: change here if to test for different inputs
    '''

    coins = ('F', 'B')
    flips = ['H', 'H', 'H', 'H', 'H', 'T', 'T', 'T', 'T', 'T']

    coin_start_probability = {
        'F': 0.5,
        'B': 0.5
    }

    # Dealer switching
    transition_probability = {
        'F': {'F': 0.9, 'B': 0.1},
        'B': {'F': 0.1, 'B': 0.9}
    }

    emission_probability = {
        'F': {'H': 0.5, 'T': 0.5},
        'B': {'H': 0.75, 'T': 0.25}
    }

    '''
    Part 1, viterbi matrix and viterbi path, leave this part as is
    '''

    print("~~PART 1:~~")
    print("Viterbi matrix:")

    # Run Viterbi
    viterbi = Viterbi(flips, coins, coin_start_probability, transition_probability, emission_probability)
    viterbi.viterbi()
    viterbi.print_viterbi_matrix()
    print()

    # Run Viterbi path
    print("Viterbi path:")
    viterbi.compute_viterbi_path()
    print()

    '''
    Part 2: backward, forward and posterior matrices
    '''

    print("~~PART 2:~~")
    print("Backward matrix:")
    viterbi.calculate_backward_probabilities()

    print("Forward matrix:")
    viterbi.calculate_forward_probabilities()

    print("Posterior:")

    # What is the probability that the "T" at the seventh position is generated by a biased coin?
    viterbi.calculate_posterior_probability(7, 'B')


if __name__ == "__main__":
    main()
