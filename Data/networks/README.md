# Networks GraphML files

This directory contains University Networks saved as GraphML files. These networks are used as input files to run python simulations.

Classroom networks are stored separately for each day of the week:<br>
M - Monday<br>
T - Tuesday<br>
W - Wednesday<br>
R - Thursday<br>
F - Friday<br>
S - Saturday<br>
U - Sunday<br>

Classroom networks without intervention are named as ClassNetwork{day}.graphml, while Classroom networks where social distancing is enforced (which results in fewer students per class and rotation of the students if the class is meeting several times a week) are stored as ClassNetPlatoon{day}.graphml files.

Student residential housing networks are stored as RoomNet.graphml, HouseholdNet.graphml, FloorNet.graphml, and Building_net.graphml files. HouseholdNet_10.graphml contains a network in which no more than 10 students per household is enforced. 
