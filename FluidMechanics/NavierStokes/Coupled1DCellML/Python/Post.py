##################
#  Post Process  #
##################

# Set the reference values
Ts = 0.001                   # Time     (s)
As = 100.0e-6                # Area     (m2)
Qs = 100.0e-6                # Flow     (m3/s)

# Set the time parameters
DYNAMIC_SOLVER_NAVIER_STOKES_STOP_TIME      = 700.0
DYNAMIC_SOLVER_NAVIER_STOKES_TIME_INCREMENT = 1.0

# Node number for data extraction
nodes = [5,9,31,87]
for node in nodes:
    # Create the Result file to store the flow and area versus time for each node
    createFile = open('Result_node'+str(node),'w+')
    # Loop through the output files
    for x in range(0,int(DYNAMIC_SOLVER_NAVIER_STOKES_STOP_TIME/DYNAMIC_SOLVER_NAVIER_STOKES_TIME_INCREMENT)):
        outputFile = open("output/MainTime_"+str(x)+".part0.exnode")
        outputLINE = outputFile.readlines()
        for i in range(0,len(outputLINE)):
            if outputLINE[i].split() == ['Node:', str(node)]:
                # Extract the variables from output files
                Flow = float(''.join(outputLINE[i+4].split()))
                Area = float(''.join(outputLINE[i+5].split()))
                A0   = float(''.join(outputLINE[i+18].split()))
                E    = float(''.join(outputLINE[i+19].split()))
                H0   = float(''.join(outputLINE[i+20].split()))
                Time = x*Ts*DYNAMIC_SOLVER_NAVIER_STOKES_TIME_INCREMENT
                Pressure = (E*H0*1.7725/A0)*(((As*Area)**0.5)-(A0**0.5))*0.0075     # Pressure(mmHg)
                Area = (Area*As-A0)*1000000.0
                Flow = Flow*Qs*1000000.0
                # Write in the Result file
                Result = open('Result_node'+str(node),'a+')
                print >> Result,"%.4f"%Time,Flow,Pressure,Area
print "."
print "."
print "."
print "Processing Completed!"
