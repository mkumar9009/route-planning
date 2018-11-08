"""Capacitated Vehicle Routing Problem with Time Windows (CVRPTW).
"""
from __future__ import print_function
from six.moves import xrange
from ortools.constraint_solver import pywrapcp
from ortools.constraint_solver import routing_enums_pb2
from math import sin, cos, sqrt, atan2, radians
import datetime
import pdb
import csv
###########################
# Problem Data Definition #
###########################
class Vehicle():
    """Stores the property of a vehicle"""
    def __init__(self):
        """Initializes the vehicle properties"""
        self._capacity = 10
        # Travel speed: 5km/h to convert in m/min
        self._speed = 25 * 60 / 3.6

    @property
    def capacity(self):
        """Gets vehicle capacity"""
        return self._capacity

    @property
    def speed(self):
        """Gets the average travel speed of a vehicle"""
        return self._speed

class ProdVol():
    """Tells the product Volume"""
    def __init__(self):
        self._cakes = 2
        self._flowers = 0.15

    @property
    def cakes(self):
        return self._cakes

    @property
    def flowers(self):
        return self._flowers

class CityBlock():
    """City block definition"""
    @property
    def width(self):
        """Gets Block size West to East"""
        return 228/2

    @property
    def height(self):
        """Gets Block size North to South"""
        return 80

class DataProblem():
    """Stores the data for the problem"""
    def __init__(self,locations,orders,order_timeslots):
        """Initializes the data for the problem"""
        self._vehicle = Vehicle()
        self._num_vehicles = 7

        # locations in meters using the city block dimension
        city_block = CityBlock()
        self._locations = [(
            (loc[0])*city_block.width,
            (loc[1])*city_block.height) for loc in locations]

        self._depot = 0
        self._demands =[]
        self._demands.insert(0,0)
        prod_vol = ProdVol()
        vehicle = Vehicle()
        for i in range(1,len(orders)):
            cakes = float(orders[i]['prod1'])
            flowers = float(orders[i]['prod2'])
            demand = cakes*prod_vol.cakes+prod_vol.flowers*flowers
            if demand > vehicle.capacity:
                print ("Order no", i, "can't be scheduled since its over capacity")
                print ("Remove Order no", i, "from list")
                exit() 
            self._demands.append(demand)

#        self._demands=[0,1,2,1,1,1,1,1,2,2,1,2]
        self._time_windows = []
        self._time_windows.append((0,0))
#            [(0, 0),(0, 8), (50, 60),  (0, 10), (10, 20), (0, 10), (75, 85),(85, 95), (5, 15), (45, 55), (30, 40),(0,0)] 
        for i in range(1,len(order_timeslots)):
            curr_time = datetime.datetime.now()
            hours_left = order_timeslots[i][1] - (curr_time.hour+1)
#            print(curr_time.hour,hours_left,order_timeslots)
            max_minutes_left   = hours_left*60 + 60-curr_time.minute
            min_time_b4_del = 60-curr_time.minute+(order_timeslots[i][0]-(curr_time.hour+1))*60
            max_time = max_minutes_left 
            if max_time < 0:
                self._time_windows.append((0,10)) #if time has passed deliver this order in first 10 minutes
            else:
                self._time_windows.append( (min_time_b4_del if min_time_b4_del>0 else 0,max_time))
        print ("Time Windows Data......................................................")
        print (self._time_windows)
#        self._time_windows = \
#            [(0, 0),(0, 8), (50, 60),  (0, 10), (10, 20), (0, 10), (75, 85),(85, 95), (5, 15), (45, 55), (30, 40)] 

    @property
    def vehicle(self):
        """Gets a vehicle"""
        return self._vehicle

    @property
    def num_vehicles(self):
        """Gets number of vehicles"""
        return self._num_vehicles

    @property
    def locations(self):
        """Gets locations"""
        return self._locations

    @property
    def num_locations(self):
        """Gets number of locations"""
        return len(self.locations)

    @property
    def depot(self):
        """Gets depot location index"""
        return self._depot

    @property
    def demands(self):
        """Gets demands at each location"""
        return self._demands

    @property
    def time_per_demand_unit(self):
        """Gets the time (in min) to load a demand"""
        return 10 # 5 minutes/unit

    @property
    def time_windows(self):
        """Gets (start time, end time) for each locations"""
        return self._time_windows

    @property
    def calculate_timeslot_window(order):
        #get the order time slot and subtract current time from the last minute of order timeslot
        curr_time = datetime.time()
        order_ts = [10,11]
        max_time = order_ts[1] - curr_time
        if max_time < 0:
            return [0,10] #if time has passed deliver this order in first 10 minutes

        return [0,max_time]
    

#######################
# Problem Constraints #
#######################
def manhattan_distance(position_1, position_2):
    """Computes the Manhattan distance between two points"""
    return (abs(position_1[0] - position_2[0]) +
            abs(position_1[1] - position_2[1]))


def geo_distance(coord_1,coord_2):
    "Computes the distance between two cordinates"
    R = 6378137 #in meters
    lat1 = radians(coord_1[0])
    lon1 = radians(coord_1[1])
    lat2 = radians(coord_2[0])
    lon2 = radians(coord_2[1])

    dlon = lon2 - lon1
    dlat = lat2 - lat1

    a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))

    return R * c    

class CreateDistanceEvaluator(object):
    """Creates callback to return distance between points."""
    def __init__(self, data):
        """Initializes the distance matrix."""
        self._distances = {}

        # precompute distance between location to have distance callback in O(1)
        for from_node in xrange(data.num_locations):
            self._distances[from_node] = {}
            for to_node in xrange(data.num_locations):
                if from_node == to_node:
                    self._distances[from_node][to_node] = 0
                else:
                    self._distances[from_node][to_node] = (
                        manhattan_distance(
                            data.locations[from_node],
                            data.locations[to_node]))

    def distance_evaluator(self, from_node, to_node):
        """Returns the manhattan distance between the two nodes"""
        return self._distances[from_node][to_node]

class CreateDemandEvaluator(object):
    """Creates callback to get demands at each location."""
    def __init__(self, data):
        """Initializes the demand array."""
        self._demands = data.demands

    def demand_evaluator(self, from_node, to_node):
        """Returns the demand of the current node"""
        del to_node
        return self._demands[from_node]

def add_capacity_constraints(routing, data, demand_evaluator):
    """Adds capacity constraint"""
    capacity = "Capacity"
    routing.AddDimension(
        demand_evaluator,
        0, # null capacity slack
        data.vehicle.capacity, # vehicle maximum capacity
        True, # start cumul to zero
        capacity)

class CreateTimeEvaluator(object):
    """Creates callback to get total times between locations."""
    @staticmethod
    def service_time(data, node):
        """Gets the service time for the specified location."""
        return data.demands[node] * data.time_per_demand_unit

    @staticmethod
    def travel_time(data, from_node, to_node):
        """Gets the travel times between two locations."""
        if from_node == to_node:
            travel_time = 0
        else:
            travel_time = manhattan_distance(
                data.locations[from_node],
                data.locations[to_node]) / data.vehicle.speed
        return travel_time

    def __init__(self, data):
        """Initializes the total time matrix."""
        self._total_time = {}
        # precompute total time to have time callback in O(1)
        for from_node in xrange(data.num_locations):
            self._total_time[from_node] = {}
            for to_node in xrange(data.num_locations):
                if from_node == to_node:
                    self._total_time[from_node][to_node] = 0
                else:
                    self._total_time[from_node][to_node] = int(
                        self.service_time(data, from_node) +
                        self.travel_time(data, from_node, to_node))

    def time_evaluator(self, from_node, to_node):
        """Returns the total time between the two nodes"""
        return self._total_time[from_node][to_node]

def add_time_window_constraints(routing, data, time_evaluator):
    """Add Global Span constraint"""
    time = "Time"
    horizon = 120
    routing.AddDimension(
        time_evaluator,
        horizon, # allow waiting time
        horizon, # maximum time per vehicle
        False, # don't force start cumul to zero since we are giving TW to start nodes
        time)
    time_dimension = routing.GetDimensionOrDie(time)
    for location_idx, time_window in enumerate(data.time_windows):
        if location_idx == 0:
            continue
        index = routing.NodeToIndex(location_idx)
        time_dimension.CumulVar(index).SetRange(time_window[0], time_window[1])
        routing.AddToAssignment(time_dimension.SlackVar(index))
    for vehicle_id in xrange(data.num_vehicles):
        index = routing.Start(vehicle_id)
        time_dimension.CumulVar(index).SetRange(data.time_windows[0][0], data.time_windows[0][1])
        routing.AddToAssignment(time_dimension.SlackVar(index))

###########
# Printer #
###########
class ConsolePrinter():
    """Print solution to console"""
    def __init__(self, data, routing, assignment):
        """Initializes the printer"""
        self._data = data
        self._routing = routing
        self._assignment = assignment

    @property
    def data(self):
        """Gets problem data"""
        return self._data

    @property
    def routing(self):
        """Gets routing model"""
        return self._routing

    @property
    def assignment(self):
        """Gets routing model"""
        return self._assignment

    def print(self):
        """Prints assignment on console"""
        # Inspect solution.
        #orders clubbed arei
        capacity_dimension = self.routing.GetDimensionOrDie('Capacity')
        time_dimension = self.routing.GetDimensionOrDie('Time')
        total_dist = 0
        total_time = 0
        for vehicle_id in xrange(self.data.num_vehicles):
            index = self.routing.Start(vehicle_id)
            plan_output = 'Route for vehicle {0}:\n'.format(vehicle_id)
            route_dist = 0
            orders_clubbed=''
            while not self.routing.IsEnd(index):
                node_index = self.routing.IndexToNode(index)
                next_node_index = self.routing.IndexToNode(
                    self.assignment.Value(self.routing.NextVar(index)))
                route_dist += manhattan_distance(
                    self.data.locations[node_index],
                    self.data.locations[next_node_index])
                load_var = capacity_dimension.CumulVar(index)
                route_load = self.assignment.Value(load_var)
                time_var = time_dimension.CumulVar(index)
                time_min = self.assignment.Min(time_var)
                time_max = self.assignment.Max(time_var)
                slack_var = time_dimension.SlackVar(index)
                slack_min = self.assignment.Min(slack_var)
                slack_max = self.assignment.Max(slack_var)
                if node_index !=0:
                    orders_clubbed = orders_clubbed + str(node_index)+','
                plan_output += ' Now Starting from order no. {0} Load({1}) Time during which it will be delivered ({2},{3}) Slack({4},{5}) \n->'.format(
                    node_index,
                    route_load,
                    time_min, time_max,
                    slack_min, slack_max)
                index = self.assignment.Value(self.routing.NextVar(index))

            node_index = self.routing.IndexToNode(index)
            load_var = capacity_dimension.CumulVar(index)
            route_load = self.assignment.Value(load_var)
            time_var = time_dimension.CumulVar(index)
            route_time = self.assignment.Value(time_var)
            time_min = self.assignment.Min(time_var)
            time_max = self.assignment.Max(time_var)
            total_dist += route_dist
            total_time += route_time
            plan_output += ' Reaching to Depot {0} Load({1}) Time ({2},{3})\n'.format(node_index, route_load, time_min, time_max)
            plan_output += 'Distance of the route: {0}m\n'.format(route_dist)
            plan_output += 'Load of the route: {0}\n'.format(route_load)
            plan_output += 'Time of the route: {0}min\n'.format(route_time)
            plan_output += 'Orders clubbed are: {0}\n'.format(orders_clubbed)
            print(plan_output)
        print('Total Distance of all routes: {0}m'.format(total_dist))
        print('Total Time of all routes: {0}min'.format(total_time))

########
# Main #
########
def main():

    with open('details.csv', mode='r') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        line_count = 0
        locations = []
        orders = []
        timeslots=[]
        for row in csv_reader:
            locations.append((float(row['lat']),float(row['long'])))
            orders.append({'prod1':row['prod1'],'prod2':row['prod2']})
            timeslots.append((int(row['start_timeslot']),int(row['end_timeslot'])))

    """Entry point of the program"""
    # Instantiate the data problem.
    data = DataProblem(locations,orders,timeslots)
    # Create Routing Model
    routing = pywrapcp.RoutingModel(data.num_locations, data.num_vehicles, data.depot)
    # Define weight of each edge
    distances = CreateDistanceEvaluator(data)
    distance_evaluator = CreateDistanceEvaluator(data).distance_evaluator
    routing.SetArcCostEvaluatorOfAllVehicles(distance_evaluator)
    # Add Capacity constraint
    demand_evaluator = CreateDemandEvaluator(data).demand_evaluator
    add_capacity_constraints(routing, data, demand_evaluator)
    # Add Time Window constraint
    time_wds= CreateTimeEvaluator(data)
    time_evaluator = CreateTimeEvaluator(data).time_evaluator
    #pdb.set_trace()
    add_time_window_constraints(routing, data, time_evaluator)

    # Setting first solution heuristic (cheapest addition).
    search_parameters = pywrapcp.RoutingModel.DefaultSearchParameters()
    search_parameters.first_solution_strategy = (
        routing_enums_pb2.FirstSolutionStrategy.PATH_CHEAPEST_ARC)
    # Solve the problem.
    assignment = routing.SolveWithParameters(search_parameters)
    printer = ConsolePrinter(data, routing, assignment)
    printer.print()

if __name__ == '__main__':
    main()
