# Dataset example

This directory contains 5 comma-separated files that describe the demographics, course registrations, schedules, and information on the student resident housing of an artificial university. The population of this university consists of 2000 undergraduate, 1000 masters and 500 PhD students as well as 250 faculty members.

These files are then converted into networks and saved as GraphML formatted files using *pre-process_classroom_data.R* and *pre-process_house_data.R* scripts located in R directory.

### course_info.csv
8 columns:<br> 
First column - courseID (character) - a unique identifier of a course;<br> 
Other seven columns  are set to TRUE if the course is taught during this day of the week and FALSE - otherwise;

### faculty_courses.csv
2 columns:<br>
courseID (character) - a unique identifier of a course;<br>
facultyID (numeric) - a unique ID for a faculty member (within university population)

### student_courses.csv
2 columns:<br>
studentID (numeric) - a unique ID for a student<br>
courseID (character) - a unique identifier of a course taken by the student

### housing_info.csv
4 columns:<br>
studentID (numeric) - a unique student ID<br>
buildingID (character) - a unique ID of a residential housing (usually building address)<br>
floor (numeric) - floor number<br>
room (numeric) - room number

### pop_info.csv
This file describes the artificial university population. In this example only students and faculty are present. <br>
7 columns:<br>
id (numeric) - student or faculty ID;<br>
age (numeric) - age; for faculty age "brackets" are used - 35 for [30-39], 45 for [40-49], etc.;<br>
sex (binary) - gender identifier (0 or 1);<br>
group (numeric) - 1 - faculty; 2 - staff; 3 - student;
campResident (numeric) - 0 - does not live on campu, 1 - lives in a small dorm, 2 - lives in a large dorm;<br>
category (numeric) - risk group ( 1 - highest risk, 4 - lowest risk). The risk group is determined by a number of factors, i.e. student is living on campus, faculty or staff meeting (or not) with students, person taking public transportation, etc. ;<br>
undergrad (binary) - 1 - undergraduate student, 0 - otherwise

Variables *group*, *campResident*, *category* and *undergrad* are used to run specific interventions.



