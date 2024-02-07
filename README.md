# Robust-Koopman-MPC


This repository contains material related  to the Automatica article "Robust tube-based model predictive control with Koopman operators"  

## Table of Contents

### Tutorials

The tutorials lead you through implementing the code files uploaded. 

* DSLC_training: Implement the code to verify the convergence condition of actor-critic learning and the closed-stability condition under actor-critic learning in the receding horizon control framework. The code is implemented in Matlab.
* DSLC_xtdrone: Deploy the control policy to control a number of multirotor drones in the Gazebo. This part is based on XTDrone, PX4, and MAVROS, containing materials related to  [XTDrone project](https://github.com/robin-shaun/XTDrone/blob/master). The code is implemented in Python. 
  * DSLC_xtdrone6: Control 6 multirotor drones to realize formation control and transformation.
  * DSLC_xtdrone18: Control 18 multirotor drones to realize formation control and transformation.
* DSLC_solving_one_robot_control: Implement the code to solve the centralized control problem of one robot distributedly and compare it with the centralized version. The code is implemented in Matlab.

## Dependencies

 1. [MPT toolbox](https://www.mpt3.org).
 2. [Yalmip toolbox](https://yalmip.github.io/).
 3. [Robust tube MPC project](https://github.com/HiroIshida/robust-tube-mpc/blob/master/example/example_tubeMPC.m), for computing the polytopic invariant sets.

## Run DSLC_xtdrone

To run the code in this repository, follow the instructions below.

1. Load worlds and drones.
    ```bash
    roslaunch multi_vehicle.launch
      ```
   
3. Obtain the position information of drones. Replace 6 with the number in the name of the selected file folder.
    ```bash
    python3 get_local_pose.py iris 6
      ```
   
5. Build the communication network among drones.
    ```bash
    multi_vehicle_communication.sh
      ```
   
6. Keyboard control code.
    ```bash
    python3 multirotor_keyboard_control_promotion.py
      ```
    *Use the keyboard to control all drones to take off and press ‘s’ to hover after a desired height. Then press ‘g’ to enter leader control mode.
   
7. Run the DSLC code for formation control.
    ```bash
    run_formation_promotion.sh
      ```
   *Note: When the script is running, and the drones are stationary, switch to the keyboard control terminal to press ‘w’ to give the leader a specified velocity. After the drones achieve the specified formation, press ‘f’ or ‘h’ to turn or press numbers 0-9 to change the formation.

8. run the baseline controller for comparison  in a straight-line formation scenario.
    ```bash
    run_formation_baseline.sh
    ```
    
9. Run the following script for the figure plot.
    ```bash
    python3 draw_figure.py
    ```
