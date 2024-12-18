import sys
# relative path to folder which contains the Sandbox module
sys.path.insert(1, '../../')
from Sandbox_V1_4 import *

from tqdm import tqdm
import networkx as nx

# import Homeostat from lab3 part 1 folder
#sys.path.insert(1, 'lab2_part1_homeostat/')
from Homeostat import *

# import local disturbances file
from homeostat_disturbances import *

def run_once(plot_data, n_units=4, dt=0.01, viability_scale=1, k=1, duration=1000, experiment=2, connection_pattern = "Fully Connected", show_snapshots = False):
    # Homeostat parameters
    upper_limit = 20
    lower_limit = -20
    upper_viability = viability_scale
    lower_viability = -viability_scale
    adapt_fun = random_val
    weights_set = None
    test_interval = 10

    # # uncomment to use a discrete weight set
    # adapt_fun = random_selector
    # weights_set = np.linspace(-1, 1, 26)

    # Simulation parameters
    t = 0
    ts = [t]
    # duration = 1000

    disturbance_start = 0

    # construct Homeostat
    homeostat = Homeostat(n_units=n_units, upper_limit=upper_limit, lower_limit=lower_limit, upper_viability=upper_viability, 
                          lower_viability=lower_viability, adapt_fun=adapt_fun, weights_set=weights_set, test_interval=test_interval, connection_pattern=connection_pattern)


    # manipulate damping parameters of springs
    for unit in homeostat.units:
        unit.k = k

    if experiment > 0:
        # start Homeostat from stable position
        for unit in homeostat.units:
            unit.thetas[-1] = 0
            unit.theta_dots[-1] = 0

    if experiment == 1:
        # create a list of disturbances
        # - positive and negative impulses alternating
        disturbances = [ImpulseDisturbanceSource(unit=homeostat.units[0], start_times=[200, 300, 400, 800], mag=10),
        ImpulseDisturbanceSource(unit=homeostat.units[0], start_times=[250, 350, 450], mag=-10)]

    elif experiment == 2:
    # create a list of disturbances
    # - a square wave disturbance which can be turned on/off
        disturbances = [SquareWaveDisturbanceSource(unit=homeostat.units[0], start_times=[50, 700], stop_times=[450], amp=10, phase_shift=10)]

    elif experiment == 3:
        # Experiment 3: Introduce disturbance at 200s and measure robustness
        disturbance_time = 200
        disturbance_magnitude = 5  # Define the magnitude of the disturbance

    else:
        disturbances = []
        
    recovery_time = duration - 200
    disturbed = False
    count = 0
    # Main Homeostat simulation loop
    while t < duration:
        homeostat.step(dt)
        if experiment == 3 and t >= disturbance_time and not disturbed:

            if show_snapshots:
                # Before Disturbance, t = 200
                plot_snapshot(homeostat, t, viability_threshold=viability_scale, file_name=f'snapshot_{int(t)}-bd.png')

            apply_disturbance(homeostat.units[0], disturbance_magnitude)
            disturbance_start = t 
            disturbed = True
            
            if show_snapshots:
                # After Disturbance, t = 200
                plot_snapshot(homeostat, t, viability_threshold=viability_scale, file_name=f'snapshot_{int(t)}-ad.png')

        t += dt
        ts.append(t)
        
        if experiment == 3 and disturbed and show_snapshots:

            if count % 2500 == 0: 
                plot_snapshot(homeostat, t, viability_threshold=viability_scale, file_name=f'snapshot_{int(t)}.png')


        if experiment == 3 and all(abs(unit.thetas[-1] - 0) < 0.1 for unit in homeostat.units) and disturbed and recovery_time == duration - 200:
            # Measure recovery time
            recovery_time = t - disturbance_start

        count += 1

    if plot_data:
        # plot Homeostat essential variables and weights over time
        fig, ax = plt.subplots(2, 1)

        # plot all homeostat unit variables over time
        for i, unit in enumerate(homeostat.units):
            ax[0].plot(ts, unit.thetas, label='Unit ' + str(i) + ': essential variable')
        ax[0].plot([ts[0], ts[-1]], [upper_viability, upper_viability], 'r--', label='upper viable boundary')
        ax[0].plot([ts[0], ts[-1]], [lower_viability, lower_viability], 'g--', label='lower viable boundary')
        ax[0].set_title('Essential variables')
        ax[0].set_xlabel('t')
        ax[0].set_ylabel('Essential variable')
        ax[0].legend()
        ax[0].legend(fontsize="6", loc ="lower right")


    if experiment == 3:
        return recovery_time
    
    return homeostat


def run_all_homeostats(n_units_array, n_runs, duration, connection_pattern, experiment=0):

    """
        Runs the homeostat n_runs for each number of units.
        :param n_units_array: Unit number array
        :param n_runs: Number of runs
        :param duration: Duration
        :param connection_pattern: Connection pattern
        :param experiment: Experiment number
    """
        
    all_homeostats = []
    dt = 0.01
        
    for i, n_units in enumerate(tqdm(n_units_array)):
        homeostat_by_run = []
        for i in range(n_runs):
            homeostat = run_once(plot_data=False, n_units=n_units, dt=dt, duration=duration, viability_scale=1, experiment=experiment, connection_pattern=connection_pattern)
            homeostat_by_run.append(homeostat)

        all_homeostats.append(homeostat_by_run)

    return (all_homeostats, connection_pattern)

def compare_all_complexities(all_pattern_homeostats, n_units_array, n_runs):

    """
        Generates graph of all complexities 
        :param all_pattern_homeostats: All homeostats, which already run
        :param n_units_array: Unit number array
        :param n_runs: Number of runs
    """

    all_average_complexities = []

    for one_pattern_homeostats in all_pattern_homeostats:
        average_complexities_one_pattern = []
        for homeostats in one_pattern_homeostats[0]:
            total_complexity = 0
            for homeostat in homeostats:
                complexity = calculate_complexity(homeostat)
                total_complexity += complexity

            average_complexities_one_pattern.append(total_complexity / n_runs)

        all_average_complexities.append((average_complexities_one_pattern, one_pattern_homeostats[1]))

    plt.figure()

    for i in range(len(all_average_complexities)):
        plt.plot(n_units_array, all_average_complexities[i][0], label = all_average_complexities[i][1])

    plt.xlabel("Number of units")
    plt.ylabel("Average Complexity")
    plt.title("Average Complexities of Connection Patterns with 10 Runs")
    plt.legend()
    plt.legend(fontsize="8", loc ="upper left")
    plt.show()


def average_time_to_stability_analysis(all_pattern_homeostats, n_units_array, n_runs, dt):

    """
        Generates graph of all adaptation times 
        :param all_pattern_homeostats: All homeostats, which already run
        :param n_units_array: Unit number array
        :param n_runs: Number of runs
        :param dt: dt
    """

    all_average_times = []

    for one_pattern_homeostats in all_pattern_homeostats:
        average_times = []
        for homeostats in one_pattern_homeostats[0]:
            adapting_times = 0
            for homeostat in homeostats:
                for i, unit in enumerate(homeostat.units):
                    adapting_times += sum(unit.testing_hist) * dt

            average_times.append(adapting_times / (n_runs*len(homeostat.units)))

        all_average_times.append((average_times, one_pattern_homeostats[1]))

    plt.figure()

    for i in range(len(all_average_times)):
        plt.plot(n_units_array, all_average_times[i][0], label = all_average_times[i][1])

    plt.xlabel("Number of units")
    plt.ylabel("Average Adapting Time")
    plt.title("Average Adapting Times of Connection Patterns with 10 Runs")
    plt.legend()
    plt.legend(fontsize="8", loc ="upper left")
    plt.show()

def robustness_analysis(connection_patterns, n_units_array, n_runs, dt, duration):

    """
        Generates graph of all recovery times 
        :param all_pattern_homeostats: All homeostats, which already run
        :param n_units_array: Unit number array
        :param n_runs: Number of runs
        :param duration: Duration
    """

    all_average_recovery_time = []
    
    for connection_pattern in tqdm(connection_patterns):
        recovery_times_one_pattern_one_unit = []
        for n_units in tqdm(n_units_array):
            total_recovery_time = 0
            for _ in range(n_runs):
                recovery_time = run_once(plot_data=False, n_units=n_units, dt=dt, duration=duration, viability_scale=1, experiment=3, connection_pattern=connection_pattern)
                total_recovery_time += recovery_time

            recovery_times_one_pattern_one_unit.append(total_recovery_time / n_runs)

        all_average_recovery_time.append((recovery_times_one_pattern_one_unit, connection_pattern))

    plt.figure()

    for i in range(len(all_average_recovery_time)):
        plt.plot(n_units_array, all_average_recovery_time[i][0], label = all_average_recovery_time[i][1])

    plt.xlabel("Number of units")
    plt.ylabel("Average Recovery Time")
    plt.title("Average Recovery Times of Connection Patterns with 10 Runs and Impulse Disturbance")
    plt.legend(fontsize="8", loc ="upper left")
    plt.show()


def calculate_complexity(homeostat):

    """
        Calculates complexity of homeostat
        :param homeostat: Homeostat
    """

    total_connections = 0
    for unit in homeostat.units:
        total_connections += len(unit.units)

    return total_connections

def apply_disturbance(unit, magnitude):
    """ 
        Applies a disturbance to the given unit. 
        :param unit: Unit to be disturbanced
        :param magnitude: Disturbance magnitude
    """
    unit.thetas[-1] = magnitude


def plot_snapshot(homeostat, t, viability_threshold=1, file_name=None):
    """ 
        Plots the current state of the Homeostat and indicates viable and non-viable units
        :param homeostat: Homeostat
        :param t: t moment
        :param viability_threshold: Viability threshold
        :param file_name: File name to save file
    """
    G = nx.Graph()
    
    for i, unit in enumerate(homeostat.units):
        viable = (unit.thetas[-1] <= viability_threshold) and (unit.thetas[-1] >= -viability_threshold)
        node_color = 'green' if viable else 'red'
        G.add_node(i, theta=unit.get_theta(), viable=viable, color=node_color)

    for i, unit in enumerate(homeostat.units):
        connected_indices = set()
        for j, connected_unit in enumerate(unit.units):
            connected_index = homeostat.units.index(connected_unit)

            if connected_index == i:
                G.add_edge(i, i)
            elif connected_index not in connected_indices:
                G.add_edge(i, connected_index, weight=unit.weights[j])
                connected_indices.add(connected_index)
    
    seed = 42 
    pos = nx.spring_layout(G, seed=seed)
    node_colors = [G.nodes[n]['color'] for n in G.nodes()]

    plt.title(f'Homeostat Snapshot at t={t:.2f}')

    nx.draw(G, pos, with_labels=True, node_color=node_colors)

    if file_name:
        plt.savefig(file_name)

    plt.subplots_adjust(top=0.9)  
    plt.show()



n_units_array = [4,5,6,7,8,9,10]
n_runs = 10
dt = 0.01
duration = 1000

robustness_analysis(["Fully Connected","Random Connectivity","Sparse Connectivity","Watts-Strogatz","Importance-based Connectivity","Barabási-Albert"], n_units_array, n_runs, dt, duration)

all_pattern_homeostats = []

all_pattern_homeostats.append(run_all_homeostats(n_units_array, n_runs, duration, "Fully Connected"))
all_pattern_homeostats.append(run_all_homeostats(n_units_array, n_runs, duration, "Random Connectivity"))
all_pattern_homeostats.append(run_all_homeostats(n_units_array, n_runs, duration, "Sparse Connectivity"))
all_pattern_homeostats.append(run_all_homeostats(n_units_array, n_runs, duration, "Watts-Strogatz"))
all_pattern_homeostats.append(run_all_homeostats(n_units_array, n_runs, duration, "Importance-based Connectivity"))
all_pattern_homeostats.append(run_all_homeostats(n_units_array, n_runs, duration, "Barabási-Albert"))

compare_all_complexities(all_pattern_homeostats,n_units_array, n_runs)
average_time_to_stability_analysis(all_pattern_homeostats,n_units_array, n_runs, dt)

plt.show()
