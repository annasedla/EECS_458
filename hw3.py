#matrix, backward, forward      posterior, main


class Viterbi:

    def __init__(self, observations, states, initial_probabilities, transition_probabilities, emission_probabilities):
        self.observations = observations
        self.states = states
        self.initial_probabilities = initial_probabilities
        self.transition_probabilities = transition_probabilities
        self.emission_probabilities = emission_probabilities
        self.viterbi_matrix = {}
        self.backward_matrix = {}
        self.forward_matrix = {}
        self.viterbi_path = [None] * len(self.observations)

    def viterbi(self):
        for state in self.initial_probabilities:
            self.viterbi_matrix[state] = [None] * (len(self.observations))
            self.viterbi_matrix[state][0] = (self.initial_probabilities[state] *
                                             self.emission_probabilities[state][self.observations[0]], "S")

        for i in xrange(len(self.observations)):
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

    def viterbi_path(self):
        previous_state = None
        probability = 0

        for state in self.states:
            if self.viterbi_matrix[state][-1][0] > probability:
                prob = self.viterbi_matrix[state][-1][0]
                previous_state = self.viterbi_matrix[state][-1][1]

        for i in xrange(len(self.observations)):
            self.viterbi_path[-(i+1)] = previous_state
            previous_state = self.viterbi_matrix[previous_state][-(i+1)][1]


def main():
    print('sup')


if __name__ == "__main__":
    main()
