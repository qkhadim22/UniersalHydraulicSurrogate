#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Author:   Qasim Khadim and Johannes Gerstmayr
# Contact:  qasim.khadim@outlook.com,qkhadim22 (Github)
# Date:     2024-06-18
# Copyright:This file is part of Exudyn. Exudyn is free software. 
# You can redistribute it and/or modify it under the terms of the Exudyn license. 
# See 'LICENSE.txt' for more details.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# %%
import numpy as np
import time

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
                    #--- DATA ACQUISITION----#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
Ts              =  1e-3              #tspan[1]-tspan[0] (Experimental data)
tEnd            =  12                #End time

createData      = False


if createData: 
    
    from Surrogate.PhysicsModels import Dynamics as PhysicsBased 
    
    Physics         = PhysicsBased(TimeStep=Ts, endTime=tEnd,    #Number of steps and time span
                                      Visualization = False, 
                                      verboseMode=1)

    start_time = time.time()
    
    inputVec        = Physics.SimulationInputVec()
 
    # Lift cylinder model
    [X1, Y1]        = Physics.LiftHydraulics(inputVec, solutionViewer = False,Plotting=True)
    LiftPhysics     = np.array([X1, Y1], dtype=object)   
    
    [X2, Y2]        = Physics.TiltHydraulics(inputVec, solutionViewer = False,Plotting=True) 
    TiltPhysics     = np.array([X2, Y2], dtype=object)

    # Save physics data 
    np.save('solution/FlexLiftPhysics', LiftPhysics)
    np.save('solution/FlexTiltPhysics', TiltPhysics)
    
    cpuTime     = time.time() - start_time
    print("--- craeting surrogate training data took: %s seconds ---" % (cpuTime))
    
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
                    #--- Standard training ----#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#X, y = datasets.load_iris(return_X_y=True)
# from sklearn.model_selection import train_test_split
# from sklearn import datasets
# from sklearn import svm

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
                    #--- CMA-ES Training ----#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#from Surrogate.Surrogate import Models

LiftSurrogate   = False
TiltSurrogate   = False


if LiftSurrogate:
    
    import cma
    from Surrogate.SurrogateModels import UniversalHydraulicSurrogate as UHS 
    
    from cmaes import CMA
    optimizer = CMA ( mean = np . zeros (12) , sigma =2)
    
    Training         = UHS(TimeStep=Ts, endTime=tEnd,    #Number of steps and time span
                                      Visualization = False, 
                                      verboseMode=1)
    
    w0_0        = np.array([1])
    w1_0        = np.array([-0.054, 1]) * 2e4  #-0.05499774560613969
    w2_0        = np.array([-1, 1]) * 1e3
    w3_0        = np.array([-1, 1]) * 5e5
    w4_0        = np.array([-1, 1]) * 10
    w5_0        = np.array([-1, 1]) * 10
    w6_0        = np.array([-1, 1]) * 10
    w7_0        = np.array([-1, 1]) * 1

    w_start     = np.concatenate([w0_0,w1_0, w2_0, w3_0,w4_0, w5_0,w6_0, w7_0])
    popsize     = len(w_start)
    
    for generation in range (100) :
        
        solutions = []
        for _ in range ( optimizer . population_size ):
            w_min       = optimizer . ask ()
            cost        = Training.ObjectiveFunc (w_min)
            solutions . append ((w_min , cost ))
            #print (f"{ generation =} { value =} {x =} ")
            
            print(f"Generation {generation} - cost: {cost:.6e}")
            
            
            
        optimizer . tell ( solutions )
        
        if optimizer . should_stop ():
            popsize = popsize * 2
            optimizer = CMA(mean=np.random.uniform(np.min(w_start), np.max(w_start), len(w_start)),
                         population_size=popsize, sigma=0.5)  
    
    
    
   

    # Initial parameters
   
    insigmax    = np.abs(w_start) / 3 
    
    
    custom_weights = np.linspace(1, 0, 3)
    custom_weights = custom_weights / np.sum(custom_weights) 
    
    options     = { #'CMA_recombination_weights': custom_weights,
                   'CMA_stds': insigmax.tolist(),
                   'verb_disp': 1,
                   'tolfun': 1e-4,
                   'maxiter': 10000}
    
    sigma       = 0.3
     
    start_time = time.time()
    result1     = cma.fmin(Training.ObjectiveFunc, w_start.tolist(), sigma, 
                           options=options)
    cpuTime     = time.time() - start_time
    print("--- Lift surrogate training took: %s seconds ---" % (cpuTime))
    
    
    xmin1       = result1[0]  
    fmin1       = result1[1]
    
    inputVec        = Training.SimulationInputVec()
    
    [X3, Y3]        = Training.SurrogateDynamics(inputVec, w_min, solutionViewer = False,Plotting=True) 

    
    LiftPara    = np.array(xmin1, dtype=object)
    np.save('solution/LiftSParaRigid', LiftPara)
    
if TiltSurrogate:
    # Initial parameters
    a_0         = np.array([1, -0.9]) * 1e3
    b0_0        = 1
    b1_0        = np.array([10, -10]) * 1
    b2_0        = np.array([20, -20]) * 1
    tau_0       = 1 / (2 * np.pi * 35)

    xstart      = np.concatenate([a_0, [b0_0], b1_0, b2_0, [tau_0]])
    insigmax    = np.abs(xstart) / 3

    options     = {'CMA_stds': insigmax.tolist()}
    sigma       = 0.5
     
    start_time = time.time()
    result2     = cma.fmin(Training.ObjectiveFunc, xstart.tolist(), sigma, options=options)
    cpuTime     = time.time() - start_time

    print("--- Tilt surrogate training took: %s seconds ---" % (cpuTime))
    
    xmin2       = result2[0]  
    fmin2       = result2[1]
    TiltPara    = np.array(xmin1, dtype=object)
    np.save('solution/TiltSParaRigid', TiltPara)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
                    #--- EVALUATION ----#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# %%
from Models.FlexibleMultibody1 import EvaluationModel


T               = 12

ns              = int(T/Ts) 
angleInit1      =   14.6        #np.deg2rad(14.6)  #Lift boom angle               
angleInit2      =  -58.8         #np.deg2rad(-58.8)     #Tilt boom angle 
LiftLoad        = 0     

Plotting        =  True

Evaluate        = EvaluationModel(nStepsTotal=ns, endTime=T,  Surrogate = False, Flexible=True,  
                                  loadFromSavedNPY=False, nModes=6, mL    = LiftLoad,verboseMode=1)

inputVec       = Evaluate.CreateInputVector( ns,  angleInit1,angleInit2 )
data           = Evaluate.ComputeModel(inputVec, solutionViewer = True) #solutionViewer: for visualization



data_array     = np.array(data, dtype=object)

if Plotting:
   Evaluate.Plotting(data_array)
   
   


#######################
#LOAD Experimetal data--

# If you change data here, then also change in ControlSignals file


###############################