# -*- coding: utf-8 -*-
"""
Created on Thu Apr 17 07:24:05 2025

@author: qkhadim22
"""

import time
import sys
sys.exudynFast = True #this variable is used to signal to load the fast exudyn module

import exudyn as exu
from exudyn.itemInterface import *
from exudyn.utilities import *
from exudyn.FEM import *
from exudyn.interactive import AnimateModes

feL                 = FEMinterface()
feT                 = FEMinterface()
SC                  = exu.SystemContainer()
mbs                 = SC.AddSystem() 
fileNameL           = 'Tet-4/LiftBoom/LiftBoom'        #File path
fileNameT           = 'Tet-4/TiltBoom/TiltBoom'        #File path


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
nModes              = 6                                    #Number of modes
loadFromSavedNPY    = True                                 #Set fasle to create a new feL for testing new mesh
LiftBoom            = False                                #Set True to visualize lift boom mode shapes
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



if not loadFromSavedNPY: 
    start_time                      = time.time()
    
    #Lift Boom
    inputFileNameStiffnessMatrixL       = 'Tet-4/LiftBoom/'+'LiftBoomStiffnessMatrix.txt'
    
    inputFileNameMassMatrixL            = 'Tet-4/LiftBoom/'+'LiftBoomMassMatrix.txt'
    inputFileNameNodalCoordinatesL      = 'Tet-4/LiftBoom/'+'LiftBoomNodes.txt'
    inputFileNameElementsL              = 'Tet-4/LiftBoom/'+'LiftBoomElements.txt'
    inputFileNameNodalMappingVectorL    = 'Tet-4/LiftBoom/'+'LiftBoomNodalMappingVector.txt' 
               
    feL.ReadStiffnessMatrixFromAnsys(inputFileNameStiffnessMatrixL, inputFileNameNodalMappingVectorL, verbose=True)
    feL.ReadMassMatrixFromAnsys(inputFileNameMassMatrixL, inputFileNameNodalMappingVectorL, verbose=True)
    feL.ReadNodalCoordinatesFromAnsys(inputFileNameNodalCoordinatesL, verbose=True)
    feL.ReadElementsFromAnsys(inputFileNameElementsL, verbose=True)                
    feL.SaveToFile(fileNameL)
    
    
    inputFileNameStiffnessMatrixT       = 'Tet-4/TiltBoom/'+'TiltBoomStiffnessMatrix.txt'
    inputFileNameMassMatrixT            = 'Tet-4/TiltBoom/'+'TiltBoomMassMatrix.txt'
    inputFileNameNodalCoordinatesT      = 'Tet-4/TiltBoom/'+'TiltBoomNode.txt'
    inputFileNameElementsT              = 'Tet-4/TiltBoom/'+'TiltBoomElement.txt'
    inputFileNameNodalMappingVectorT    = 'Tet-4/TiltBoom/'+'TiltBoomNodalMappingVector.txt'     
    
    feT.ReadStiffnessMatrixFromAnsys(inputFileNameStiffnessMatrixT, inputFileNameNodalMappingVectorT, verbose=True)
    feT.ReadMassMatrixFromAnsys(inputFileNameMassMatrixT, inputFileNameNodalMappingVectorT, verbose=True)
    feT.ReadNodalCoordinatesFromAnsys(inputFileNameNodalCoordinatesT, verbose=True)
    feT.ReadElementsFromAnsys(inputFileNameElementsT, verbose=True)                
    feT.SaveToFile(fileNameT)
  
else:       
    print('importing ANSYS FEM data structure of Lift Boom...')
    start_time = time.time()
    feL.LoadFromFile(fileNameL)
    feT.LoadFromFile(fileNameT)
    cpuTime = time.time() - start_time
    print("--- importing FEM data took: %s seconds ---" % (cpuTime))
    

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  

if LiftBoom:  
    print("Compute lift boom free modes... ")
    feL.ComputeEigenmodes(nModes, excludeRigidBodyModes = 6, useSparseSolver = True)
    print("Lift Boom eigen freq.=", feL.GetEigenFrequenciesHz())  

    cms1     = ObjectFFRFreducedOrderInterface(feL)   
    objFFRF1 = cms1.AddObjectFFRFreducedOrder(mbs, positionRef=[0,0,0],initialVelocity=[0,0,0],initialAngularVelocity=[0,0,0],
                                                gravity=[0,-9.81,0],color=[0.1,0.9,0.1,1.],)             
    
    mbs.Assemble()
        
    SC.visualizationSettings.nodes.show = False
    SC.visualizationSettings.openGL.showFaceEdges = True
    SC.visualizationSettings.openGL.multiSampling=4   
    SC.visualizationSettings.window.renderWindowSize = [1600,1080]
    SC.visualizationSettings.contour.outputVariable = exu.OutputVariableType.DisplacementLocal
    SC.visualizationSettings.contour.outputVariableComponent = 0 #component
            
    SC.visualizationSettings.general.autoFitScene = False #otherwise, model may be difficult to be moved
    nodeNumber = objFFRF1['nGenericODE2'] #this is the node with the generalized coordinates
    AnimateModes(SC, mbs, nodeNumber)
    import sys
    sys.exit()  
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
else:
    print("Compute tilt boom free modes... ")
    feT.ComputeEigenmodes(nModes, excludeRigidBodyModes = 6, useSparseSolver = True)
    print("Tilt Boom eigen freq.=", feT.GetEigenFrequenciesHz())  
    
    cms2     = ObjectFFRFreducedOrderInterface(feT)   
    objFFRF2 = cms2.AddObjectFFRFreducedOrder(mbs, positionRef=[0,0,0],initialVelocity=[0,0,0],initialAngularVelocity=[0,0,0],
                                                    gravity=[0,-9.81,0],color=[0.1,0.9,0.1,1.],)             
    
    mbs.Assemble()
        
    SC.visualizationSettings.nodes.show = False
    SC.visualizationSettings.openGL.showFaceEdges = True
    SC.visualizationSettings.openGL.multiSampling=4   
    SC.visualizationSettings.window.renderWindowSize = [1600,1080]
    SC.visualizationSettings.contour.outputVariable = exu.OutputVariableType.DisplacementLocal
    SC.visualizationSettings.contour.outputVariableComponent = 0 #component
        
    SC.visualizationSettings.general.autoFitScene = False #otherwise, model may be difficult to be moved
    nodeNumber = objFFRF2['nGenericODE2'] #this is the node with the generalized coordinates
    AnimateModes(SC, mbs, nodeNumber)
    import sys
    sys.exit()  