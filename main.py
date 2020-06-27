import math
import numpy as np
from matplotlib import pyplot as plt


#all distances are in micrometres
#all time units are seconds

#begin input parameters
plane_width = 30 #receiver width
plane_height = 30 #receiver height
plane_width_transmitter = 30 #transmitter width
plane_height_transmitter = 30 #transmitter height
distance_between_planes = 20
trigger_threshold = 3000 #trigger threshold value to emit molecule
time_to_relax = 1
timestamp = pow(10,-4)
D = 79.4
standart_deviation = math.sqrt(2*D*timestamp)
receiver_diameter = 5
transmitter_diameter = 5
receivers_start_point = [distance_between_planes,0,0]
transmitters_start_point = [0,0,0]
number_of_molecule = 1000 #number of molecules emitted from each cell
#end input parameters

number_of_receiver = int(plane_height/(2*receiver_diameter))*int(plane_width/(2*receiver_diameter))
number_of_transmitter = int(plane_height_transmitter/(2*transmitter_diameter))*int(plane_width_transmitter/(2*transmitter_diameter))

receivers = np.zeros((number_of_receiver,5))
transmitters = np.zeros((number_of_transmitter,5))
a_molecules = np.zeros((1,3))
a_molecules = np.delete(a_molecules, (0), axis=0)
b_molecules = np.zeros((1,3))
b_molecules = np.delete(b_molecules, (0), axis=0)

z_coordinates = []
z = receivers_start_point[2] + receiver_diameter
z_of_end_point = receivers_start_point[2] + plane_width
while (z < z_of_end_point):
    z_coordinates.append(z)
    z = z + 2 * receiver_diameter

y_coordinates = []
y = receivers_start_point[1] + receiver_diameter
y_of_end_point = receivers_start_point[1] + plane_height
while (y < y_of_end_point):
    y_coordinates.append(y)
    y = y + 2 * receiver_diameter

i = 0
for cor_y in y_coordinates:
    for cor_z in z_coordinates:
        receivers[i] = [receivers_start_point[0], cor_y, cor_z, 0, 0]
        i = i + 1


z_coordinates = []
z = transmitters_start_point[2] + transmitter_diameter
z_of_end_point = transmitters_start_point[2] + plane_width_transmitter
while(z < z_of_end_point):
    z_coordinates.append(z)
    z = z + 2*transmitter_diameter

y_coordinates = []
y = transmitters_start_point[1] + transmitter_diameter
y_of_end_point = transmitters_start_point[1] + plane_height_transmitter
while (y < y_of_end_point):
    y_coordinates.append(y)
    y = y + 2 * transmitter_diameter

i = 0
for cor_y in y_coordinates:
    for cor_z in z_coordinates:
        transmitters[i]= [transmitters_start_point[0], cor_y, cor_z, 0, 0]
        i = i+1
transmitters[0] = [transmitters_start_point[0], y_coordinates[0], z_coordinates[0], 0,0]

for tran in transmitters:
    new_molecules = np.zeros((number_of_molecule, 3), dtype=int)
    new_molecules[:, 0] = tran[0] + transmitter_diameter
    new_molecules[:, 1] = tran[1]
    new_molecules[:, 2] = tran[2]
    a_molecules = np.vstack((a_molecules, new_molecules))


#simulation begin

receivedA = []
receivedB = []
curr_received_a = 0
cum_received_a = 0
curr_received_b = 0
cum_received_b = 0
sqrt_transmitter_radius = transmitter_diameter ** 2
sqrt_receiver_radius = receiver_diameter ** 2


iter = 1
ind_col = 0
for i in range(1, int(time_to_relax / timestamp + 1)):
    len_a_molecules = len(a_molecules)
    len_b_molecules = len(b_molecules)
    is_a_zero = False
    if len_a_molecules == 0:
        is_a_zero = True
    is_b_zero = False
    if len_b_molecules == 0:
        is_b_zero = True


    dist_a =  np.random.normal(0, standart_deviation, (len_a_molecules, 3))+ a_molecules
    filter_validation_a_molecules = np.array(np.ones(len_a_molecules), dtype=bool)
    filter_received_a_molecules = np.array(np.zeros(len_a_molecules), dtype=bool)

    dist_b = np.random.normal(0, standart_deviation, (len_b_molecules, 3)) + b_molecules
    filter_validation_b_molecules = np.array(np.ones(len_b_molecules), dtype=bool)
    filter_received_b_molecules = np.array(np.zeros(len_b_molecules), dtype=bool)
    #validation for a
    #received for b

    for t in transmitters:
        x = -t[0]
        y = -t[1]
        z = -t[2]
        if is_a_zero == False:
            #a molecules
            curr_molecules_a = np.column_stack((dist_a, np.zeros((len_a_molecules))))
            curr_molecules_a[:, 0] += x
            curr_molecules_a[:, 1] += y
            curr_molecules_a[:, 2] += z
            curr_molecules_a = np.power(curr_molecules_a, 2)
            curr_molecules_a[:, 3] = curr_molecules_a[:, 0] + curr_molecules_a[:, 1] + curr_molecules_a[:, 2]
            index_of_not_move = np.where(curr_molecules_a[:, 3] < sqrt_transmitter_radius)[0]
            for index in index_of_not_move:
                filter_validation_a_molecules[index] = False

        if is_b_zero == False:
            #b molecules
            curr_molecules_b = np.column_stack((dist_b, np.zeros((len_b_molecules))))
            curr_molecules_b[:, 0] += x
            curr_molecules_b[:, 1] += y
            curr_molecules_b[:, 2] += z
            curr_molecules_b = np.power(curr_molecules_b, 2)
            curr_molecules_b[:, 3] = curr_molecules_b[:, 0] + curr_molecules_b[:, 1] + curr_molecules_b[:, 2]
            index_of_received = np.where(curr_molecules_b[:, 3] < sqrt_transmitter_radius)[0]
            for index in index_of_received:
                filter_received_b_molecules[index] = True
            t[4] = t[4] + len(index_of_received)
            curr_received_b += len(index_of_received)
            cum_received_b += len(index_of_received)



    for r in receivers:
        x = -r[0]
        y = -r[1]
        z = -r[2]
        if is_a_zero == False:
            #a
            curr_molecules_a = np.column_stack((dist_a, np.zeros((len_a_molecules))))
            curr_molecules_a[:, 0] += x
            curr_molecules_a[:, 1] += y
            curr_molecules_a[:, 2] += z
            curr_molecules_a = np.power(curr_molecules_a, 2)
            curr_molecules_a[:, 3] = curr_molecules_a[:, 0] + curr_molecules_a[:, 1] + curr_molecules_a[:, 2]
            index_of_received = np.where(curr_molecules_a[:, 3] < sqrt_receiver_radius)[0]
            for index in index_of_received:
                filter_received_a_molecules[index] = True
            r[4] = r[4] + len(index_of_received)
            curr_received_a += len(index_of_received)
            cum_received_a += len(index_of_received)

        if is_b_zero == False:
            #b
            curr_molecules_b = np.column_stack((dist_b, np.zeros((len_b_molecules))))
            curr_molecules_b[:, 0] += x
            curr_molecules_b[:, 1] += y
            curr_molecules_b[:, 2] += z
            curr_molecules_b = np.power(curr_molecules_b, 2)
            curr_molecules_b[:, 3] = curr_molecules_b[:, 0] + curr_molecules_b[:, 1] + curr_molecules_b[:, 2]
            index_of_not_move = np.where(curr_molecules_b[:, 3] < sqrt_receiver_radius)[0]
            for index in index_of_not_move:
                filter_validation_b_molecules[index] = False


    if is_a_zero == False:
        filter_valid_not_received_a = filter_validation_a_molecules & np.invert(filter_received_a_molecules)
        not_valid_a_molecules = a_molecules[np.invert(filter_validation_a_molecules)]
        not_received_a_molecules = dist_a[filter_valid_not_received_a]
        a_molecules = np.vstack((not_valid_a_molecules, not_received_a_molecules))
    receivedA.append(curr_received_a)

    if is_b_zero == False:
        filter_valid_not_received_b = filter_validation_b_molecules & np.invert(filter_received_b_molecules)
        not_valid_b_molecules = b_molecules[np.invert(filter_validation_b_molecules)]
        not_received_b_molecules = dist_b[filter_valid_not_received_b]
        b_molecules = np.vstack((not_valid_b_molecules, not_received_b_molecules))
    receivedB.append(curr_received_b)

    curr_received_a = 0
    curr_received_b = 0

    if cum_received_a > trigger_threshold:
        for rec in receivers:
            new_molecules = np.zeros((number_of_molecule, 3), dtype=int)
            new_molecules[:, 0] = rec[0] - receiver_diameter
            new_molecules[:, 1] = rec[1]
            new_molecules[:, 2] = rec[2]
            b_molecules = np.vstack((b_molecules, new_molecules))
        cum_received_a = 0
    if cum_received_b > trigger_threshold:
        for tran in transmitters:
            new_molecules = np.zeros((number_of_molecule, 3), dtype=int)
            new_molecules[:, 0] = tran[0] + transmitter_diameter
            new_molecules[:, 1] = tran[1]
            new_molecules[:, 2] = tran[2]
            a_molecules = np.vstack((a_molecules, new_molecules))
        cum_received_b = 0

    ind_col += 1

fig, (ax1) = plt.subplots(1)

merge = (time_to_relax/timestamp)/100

newReceivedA = np.zeros(merge+ 1)
for i in range(merge):
    subList = receivedA[i*100:(i+1)*100]
    newReceivedA[i+1] = sum(subList)
newReceivedB = np.zeros(merge+1)
for i in range(merge):
    subList = receivedB[i*100:(i+1)*100]
    newReceivedB[i+1] = sum(subList)

ax1.plot(np.arange(merge+1), newReceivedA)
ax1.plot(np.arange(merge+1), newReceivedB)
ax1.legend(["a", "b"])
plt.xlabel('time(s)*10^-2')
plt.ylabel('Number of received molecules per timestamp')
plt.show()

