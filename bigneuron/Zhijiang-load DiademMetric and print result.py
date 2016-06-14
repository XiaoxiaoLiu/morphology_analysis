'''
Created on 2016/06/11

@author: Zhijiang-PC
'''

import jpype
import os.path
 
classpath ='C:/vaa3d/DiademMetric.jar'

jpype.startJVM(jpype.getDefaultJVMPath(),"-ea","-Djava.class.path=%s"%(classpath))
ReadSWC_class = jpype.JClass("org.krasnow.cng.data.ReadSWC")
SwcTreeNode_class=jpype.JClass("org.krasnow.cng.domain.SwcTreeNode")

DiademMetric_Class=jpype.JClass("org.krasnow.cng.diadem.DiademMetric")

File_Class= jpype.JClass("java.io.File")
List_Class=jpype.java.util.ArrayList()
t = ReadSWC_class()     
Swc_object=SwcTreeNode_class()

a_object=File_Class("C:/vaa3d/DiademMetric/ExampleGoldStandard.swc")
b_object=File_Class("C:/vaa3d/DiademMetric/ExampleTest.swc")
DiademMetric_object=DiademMetric_Class(a_object,b_object)
a_value =ReadSWC_class.convertSwcToBinaryTreeList(a_object).size();
b_value=DiademMetric_object.getXYThreshold()

jprint = jpype.java.lang.System.out.println
jprint(a_value)
jprint(b_value)
jpype.shutdownJVM()