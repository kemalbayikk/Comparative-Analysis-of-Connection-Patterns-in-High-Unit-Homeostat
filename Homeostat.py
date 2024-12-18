import sys
import networkx as nx
# relative path to folder which contains the Sandbox module
sys.path.insert(1, '../../')
from Sandbox_V1_4 import *

import copy as cp

'''
    A class to simulate a single unit in a simulation of Ashby's Homeostat machine.
'''
class Unit(System):

    def __init__(self, test_interval, adapt_fun, upper_viability=1, lower_viability=-1, upper_limit=np.Inf, lower_limit=-np.Inf, m=1, k=1, l=1, p=2, q=1, theta0=2, theta_dot0=0, weights_set=None, adapt_enabled=True):

        # Ashby's homeostat used a discrete set of weights. If we want
        # to do that, we need to pass a list of weights into here, and
        # our adapt_fun needs to be written to work with them (e.g. see
        # random_selector)
        self.weights_set = weights_set
        # minimum time to wait between changing the weights
        self.test_interval = test_interval
        # clock
        self.t = 0
        # system equation parameters
        self.m = m # moment of inertia
        self.k = k # damping due to friction
        self.l = l # valve gain
        self.p = p # electric potential at one end of trough
        self.q = q # electric potential at the other end of trough
        # "physical" upper limit for theta - the variable (or needle) can't move past this point
        self.upper_limit = upper_limit
        # "physical" lower limit for theta - the variable (or needle) can't move past this point
        self.lower_limit = lower_limit
        # upper limit for viability
        self.upper_viability = upper_viability
        # lower limit for viability
        self.lower_viability = lower_viability

        self.thetas = [theta0] # the system variable. theta is used instead of x, because in SituSim x is always used to mean a spatial coordinate
        self.theta_dots = [theta_dot0] # the system state is [theta, theta_dot]
        self.theta_dotdots = [0] # the input to the system is converted into an acceleration (technically, it is a force, but mass is not modelled here)
        self.units = [] # a list of connected Units
        self.weights = [] # a list of weights, which are applied to connected units
        self.weights_hist = [] # a record of weights over time
        self.testing = False # this boolean variable keeps track of whether the system is in the process of adaptation
        self.testing_hist = [self.testing] # a record of when the unit is adapting
        self.inputs_hist = [] # a record of all inputs to a unit over time, including the feedback from the unit itself
        self.adapt_fun = adapt_fun #
        self.self_ind = 0 # used to be able to set the unit's feedback connection to always be negative
        self.test_times = [] # times when the unit starts testing new weights
        self.timer = 0 # a timer for how long to wait after changing a unit's weights
        self.adapt_enabled = adapt_enabled # Unit will only adapt if this variable is set to True

    # this method should be called after a unit's connections have
    # been made, and before it is run
    def initialise(self):
        self.weights_hist.append(self.weights)

    # randomise the parameters which affect the unit's dynamics
    # (i.e. all parameters *apart from* the connection weights)
    def randomise_params(self):

        self.m = random_in_interval(minimum=0.1, maximum=2)
        self.k = random_in_interval(minimum=0.1, maximum=2)
        self.l = random_in_interval(minimum=0.1, maximum=2)
        self.q = random_in_interval(minimum=0.1, maximum=1)
        self.p = self.q + random_in_interval(minimum=0.1, maximum=2)

        print("***** Randomising homeostat unit params *****")
        print("m", self.m)
        print("k", self.k)
        print("l", self.l)
        print("p", self.p)
        print("q", self.q)

    # step unit forwards in time
    def step(self, dt):
        # get weighted sum of inputs
        input_sum = self.update_inputs(dt)
        # integrate the system's dynamics
        self.integrate(dt, input_sum)

        # manage timer
        if self.timer > self.test_interval:
            # reset timer
            self.timer = 0
            self.testing = False

        # if not viable, then attempt to adapt
        if self.adapt_enabled and not self.test_viability():
            # not viable
            # if already testing new parameters, then continue to test
            # otherwise, try some new weights
            if not self.testing:
                self.adjust_weights(dt)
                # start testing
                self.testing = True
                # keep record of times when new weights are set
                self.test_times.append(self.t)
                # reset timer
                self.timer = 0

        # keep record of when the system is testing new weights
        self.testing_hist.append(self.testing)
        # keep history of weights over time
        self.weights_hist.append(self.weights)

        # increment clock
        self.t += dt
        # increment test timer
        self.timer += dt

        # return current state
        return self.thetas[-1]

    # get inputs from all connect units, including this unit
    def update_inputs(self, dt):

        # calculate weighted sum of inputs from all connected Units (including feedback from this Unit)
        input_sum = 0
        inputs = []
        for unit, weight in zip(self.units, self.weights):
            inputs.append(unit.get_theta() * weight)
            input_sum += inputs[-1]
        # keep record of inputs
        self.inputs_hist.append(inputs)

        return input_sum

    # integrate system's dynamics
    def integrate(self, dt, input_sum):

        # integrate the system, from acceleration to position
        theta_dotdot = (-self.k * self.theta_dots[-1]) + (self.l * (self.p - self.q) * input_sum) # calculate acceleration
        # - we integrate twice here, because this is a second order system
        theta_dot = (self.theta_dots[-1] + (self.theta_dotdots[-1] * dt)) # integrate acceleration to get velocity
        theta = (self.thetas[-1] + (self.theta_dots[-1] * dt)) # integrate velocity to get position

        # in Ashby's Homeostat, there were hard limits to how far the needle
        # (system variable) could move in either direction - enforce these limits
        # If you don't want to enforce limits, you can leave them at them
        # at the default values of +-np.Inf
        if theta > self.upper_limit:
            theta = self.upper_limit
            theta_dot = 0
        elif theta < self.lower_limit:
            theta = self.lower_limit
            theta_dot = 0

        # store system variable and its first and second derivatives
        self.thetas.append(theta)
        self.theta_dots.append(theta_dot)
        self.theta_dotdots.append(theta_dotdot)

    # test to see whether the system's essential variable (theta) is
    # within the chosen limits for viability
    def test_viability(self):

        # test whether  the Unit is within viable limits
        if ((self.thetas[-1] > self.upper_viability) or
            (self.thetas[-1] < self.lower_viability)):
            return False # return False for not viable
        return True # return True for viable

    # adjust weights
    def adjust_weights(self, dt):

        # adjust parameters
        self.weights = self.adapt_fun(dt, self.inputs_hist, self.weights_hist, self.thetas, self.theta_dots, weights_set=self.weights_set, self_ind=self.self_ind)

    # connect a Unit to this Unit
    def add_connection(self, unit, weight):
        # add Unit to list
        self.units.append(unit)
        # add connection weight to weights list
        self.weights.append(weight)

        # if adding a connection from the unit to itself, keep
        # track of the index of the weight
        if id(unit) == id(self):
            self.self_ind = len(self.units) - 1

    # get the state of the Unit variable (the full state of a Unit is
    # actually [theta, theta_dot], but other Units can only "see" theta)
    def get_theta(self):
        return self.thetas[-1]

'''
    A class to simulate Ashby's Homeostat machine.
'''
class Homeostat(System):

    def __init__(self, n_units, upper_viability, lower_viability, adapt_fun, upper_limit=np.Inf, 
                 lower_limit=-np.Inf, weights_set=None, test_interval=10, adapt_enabled=True, connection_pattern = "Sparse Connectivity"):

        # set up units
        self.units = []
        for _ in range(n_units):
            self.units.append(Unit(test_interval=test_interval, adapt_fun=adapt_fun, upper_limit=upper_limit, lower_limit=lower_limit, upper_viability=upper_viability, lower_viability=lower_viability, weights_set=weights_set, adapt_enabled=adapt_enabled))

        if connection_pattern == "Fully Connected":
            for unit in self.units:
                for unit2 in self.units:
                    if unit != unit2:
                        unit.add_connection(unit2, random_in_interval(minimum=-1, maximum=1))

            for i, unit in enumerate(self.units):
                self_weight = -np.abs(np.random.uniform(0.5, 1.5))
                unit.add_connection(unit, self_weight)


        elif connection_pattern == "Sparse Connectivity":
            num_connections = int(n_units * 0.3)
            for unit in self.units:
                possible_targets = [other_unit for other_unit in self.units if other_unit != unit]
                connected_units = random.sample(possible_targets, min(num_connections, len(possible_targets)))
                for target in connected_units:
                    weight = random.uniform(-1, 1)
                    unit.add_connection(target, weight)

            for i, unit in enumerate(self.units):
                self_weight = -np.abs(np.random.uniform(0.5, 1.5))
                unit.add_connection(unit, self_weight)

        elif connection_pattern == "Random Connectivity":
            
            connection_matrix = np.random.rand(n_units, n_units) < 0.5
            np.fill_diagonal(connection_matrix, False)
            n_units = len(self.units)
            for i in range(n_units):
                for j in range(n_units):
                    if connection_matrix[i, j]:
                        weight = np.random.uniform(-1, 1)
                        self.units[i].add_connection(self.units[j], weight)

            for i, unit in enumerate(self.units):
                self_weight = -np.abs(np.random.uniform(0.5, 1.5))
                unit.add_connection(unit, self_weight)

        elif connection_pattern == "Barabási-Albert":
            initial_connections = 2
            ba_graph = nx.barabasi_albert_graph(n_units, initial_connections)
    
            for i, unit in enumerate(self.units):
                for j in ba_graph.neighbors(i):
                    if i != j:
                        weight = np.random.uniform(-1, 1)
                        unit.add_connection(self.units[j], weight)
                
                self_weight = -np.abs(np.random.uniform(0.5, 1.5))
                unit.add_connection(unit, self_weight)

        elif connection_pattern == "Importance-based Connectivity":
            importance_scores = np.random.rand(n_units)
            total_importance = np.sum(importance_scores)
            connection_probabilities = importance_scores / total_importance

            num_connections = np.random.poisson(lam=2, size=n_units)
            num_connections = np.clip(num_connections, 1, n_units-1)

            for i, unit in enumerate(self.units):
                possible_targets = [other_unit for other_unit in self.units if other_unit != unit]
                normalized_probabilities = [connection_probabilities[j] for j in range(n_units) if self.units[j] != unit]

                normalized_probabilities /= np.sum(normalized_probabilities)

                num_conns = num_connections[i]
                weights = np.random.uniform(-1, 1, size=num_conns)
                connections = np.random.choice(possible_targets, size=num_conns, replace=False, p=normalized_probabilities)
                
                for target, weight in zip(connections, weights):
                    unit.add_connection(target, weight)

                self_weight = -np.abs(np.random.uniform(0.5, 1.5))
                unit.add_connection(unit, self_weight)
            
        elif connection_pattern == "Watts-Strogatz":

            neighborhood = 2
            rewire_prob = 0.3
            sw_graph = nx.watts_strogatz_graph(n_units, neighborhood, rewire_prob)
        
            for i, unit in enumerate(self.units):
                for j in sw_graph.neighbors(i):
                    if i != j:
                        weight = np.random.uniform(-1, 1)
                        unit.add_connection(self.units[j], weight)
                
                self_weight = -np.abs(np.random.uniform(0.5, 1.5))
                unit.add_connection(unit, self_weight)


        # initialise units
        self.initialise()
        # initialise Homeostat's time variable
        self.t = 0

    # step unit forwards in time
    def step(self, dt):
        for unit in self.units:
            unit.step(dt)

        self.t += dt

    # initialise units
    def initialise(self):
        for unit in self.units:
            unit.initialise()

    # randomise unit parameters (this does not modify connection weights)
    def randomise_params(self):
        for unit in self.units:
            unit.randomise_params()

    def draw_network(self):
        G = nx.DiGraph()
        for i, unit in enumerate(self.units):
            G.add_node(i, theta=unit.get_theta())
            for target, weight in zip(unit.units, unit.weights):
                G.add_edge(self.units.index(unit), self.units.index(target), weight=weight)

        pos = nx.spring_layout(G)
        nx.draw(G, pos, with_labels=True, node_color='skyblue', node_size=300, edge_color='k', connectionstyle='arc3,rad=0.1')
        edge_labels = nx.get_edge_attributes(G, 'weight')
        nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)
        plt.title("Homeostat Network Visualization")
        plt.show()

'''
	An example adaptation function, loosely based on Ashby's random step change.
	It differs from Ahsby's mechanism in that it does
	not choose weights from a discrete set of values, but randomly
	draws values from a uniform interval.
'''
def random_val(dt, inputs_hist, weights_hist, thetas, theta_dots, weights_set=[], self_ind = None):
    weights = []
    for _ in range(len(weights_hist[0])):
        weights.append(random_in_interval(-1, 1))

    # this ensures that the self-connection on a Unit will have a
	# negative weight - this is not a requirement, but will lead to
	# stability being found quicker
    if self_ind is not None:
        weights[self_ind] = - np.abs(weights[self_ind])

    # return new weights
    return weights

'''
	An example adaptation function, which moves weights by a small
	random amount from their current values.
'''
def random_creeper(dt, inputs_hist, weights_hist, thetas, theta_dots, weights_set=[], self_ind = None):

    weights = []
    for i in range(len(weights_hist[-1])):
        weights.append(weights_hist[-1][i] + random_in_interval(-0.05, 0.05))

    # this ensures that the self-connection on a Unit will have a
	# negative weight - this is not a requirement, but will lead to
	# stability being found quicker
    if self_ind:
        weights[self_ind] = - np.abs(weights[self_ind])

    # return new weights
    return weights

'''
	An example adaptation function, loosely based on Ashby's random step change.
	Like Ashby's mechanism it chooses weights from a discrete set of values, "weights_set".
'''
def random_selector(dt, inputs_hist, weights_hist, thetas, theta_dots, weights_set, self_ind = None):

    assert weights_set is not None, "To use random_selector, you must pass a list of weights to __init__ when you construct your Homeostat"

    weights = []
    for _ in range(len(weights_hist[0])):
        weights.append(np.random.choice(weights_set))

    # this ensures that the self-connection on a Unit will have a
	# negative weight - this is not a requirement, but will lead to
	# stability being found quicker
    if self_ind:
        weights[self_ind] = - np.abs(weights[self_ind])

    # return new weights
    return weights
